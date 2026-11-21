/*
 * Copyright 2026 Tom van Dijk
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * End-to-end test of the interactive screens of bench_solvers on a pseudo
 * terminal (POSIX only): it starts bench_solvers on a pty, presses keys and
 * clicks, and checks the screen, the exit status and that the terminal is
 * restored. The screens themselves are unit tested without a terminal; this
 * test covers the terminal I/O around them.
 *
 * Usage: test_pty <bench_solvers> <oink> <directory of test games> <work directory>
 */

#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <fcntl.h>
#include <poll.h>
#include <signal.h>
#include <sys/ioctl.h>
#include <sys/wait.h>
#include <termios.h>
#include <unistd.h>

#include "common/terminal.hpp"

namespace fs = std::filesystem;

static int failures = 0;

static void
check(const std::string& name, bool ok)
{
    if (!ok) {
        std::cerr << "FAIL [" << name << "]" << std::endl;
        failures++;
    }
}

/**
 * A program running on a pseudo terminal of 100 by 30.
 *
 * The program does not lead the session of the terminal: a keeper process does, runs the
 * program as its child, reports its exit status and keeps the terminal until the test has
 * read all output. When a session leader exits, the system revokes its terminal, and on FreeBSD
 * output that was not read yet, such as the last lines before the exit, can then be lost.
 */
class Session
{
public:
    Session(const std::string& program, const std::vector<std::string>& args, const std::string& directory)
    {
        master = posix_openpt(O_RDWR | O_NOCTTY);
        if (master < 0 || grantpt(master) != 0 || unlockpt(master) != 0) throw std::runtime_error("cannot open a pseudo terminal");
        slave_name = ptsname(master);
        winsize size{};
        size.ws_col = 100;
        size.ws_row = 30;
        ioctl(master, TIOCSWINSZ, &size);
        // keep the slave open, to read its settings after the program ends
        slave = open(slave_name.c_str(), O_RDWR | O_NOCTTY);
        tcgetattr(slave, &before);

        // the keeper reports the exit code of the program on <status>, and exits when the test
        // writes a byte to <release>; the programs of later sessions must not inherit them
        int status[2], release[2];
        if (pipe(status) != 0 || pipe(release) != 0) throw std::runtime_error("cannot create pipes");
        for (const int fd : {status[0], status[1], release[0], release[1]}) fcntl(fd, F_SETFD, FD_CLOEXEC);

        pid = fork();
        if (pid == 0) {
            close(status[0]);
            close(release[1]);
            setsid();
            const int tty = open(slave_name.c_str(), O_RDWR);
            ioctl(tty, TIOCSCTTY, 0);
            dup2(tty, 0);
            dup2(tty, 1);
            dup2(tty, 2);
            int code = 126;
            if (chdir(directory.c_str()) == 0) {
                // the signals of the keyboard are for the program
                for (const int sig : {SIGINT, SIGQUIT, SIGTSTP}) signal(sig, SIG_IGN);
                const pid_t child = fork();
                if (child == 0) {
                    for (const int sig : {SIGINT, SIGQUIT, SIGTSTP}) signal(sig, SIG_DFL);
                    setenv("TERM", "xterm-256color", 1);
                    unsetenv("NO_COLOR");
                    std::vector<char*> argv;
                    argv.push_back(const_cast<char*>(program.c_str()));
                    for (const auto& arg : args) argv.push_back(const_cast<char*>(arg.c_str()));
                    argv.push_back(nullptr);
                    execv(program.c_str(), argv.data());
                    _exit(127);
                }
                code = 127;
                int wstatus = 0;
                if (child > 0) {
                    while (waitpid(child, &wstatus, 0) < 0 && errno == EINTR) {}
                    code = WIFEXITED(wstatus) ? WEXITSTATUS(wstatus) : -1;
                }
            }
            if (write(status[1], &code, sizeof(code)) != (ssize_t)sizeof(code)) _exit(1);
            char byte;
            while (read(release[0], &byte, 1) < 0 && errno == EINTR) {}
            _exit(0);
        }
        close(status[1]);
        close(release[0]);
        status_fd = status[0];
        release_fd = release[1];
    }

    ~Session()
    {
        if (pid > 0 && !exited) {
            // the process group of the keeper: the keeper and the program
            if (kill(-pid, SIGKILL) != 0) kill(pid, SIGKILL);
            waitpid(pid, nullptr, 0);
        }
        close(status_fd);
        close(release_fd);
        close(slave);
        close(master);
    }

    /** Read what the program wrote, for at most <seconds>. */
    void read_for(double seconds)
    {
        const auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
        while (std::chrono::steady_clock::now() < end) {
            pollfd fd{master, POLLIN, 0};
            if (poll(&fd, 1, 20) <= 0) continue;
            char data[4096];
            const ssize_t n = read(master, data, sizeof(data));
            if (n <= 0) return; // the terminal is closed
            output.append(data, (size_t)n);
        }
    }

    /** Wait until the screen shows <text>; false after a timeout. */
    bool expect(const std::string& text, double seconds = 20)
    {
        const auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
        while (std::chrono::steady_clock::now() < end) {
            if (ui::strip(output.substr(seen)).find(text) != std::string::npos) {
                seen = output.size();
                return true;
            }
            read_for(0.05);
        }
        std::cerr << "timeout waiting for \"" << text << "\"; screen:\n" << ui::strip(output.substr(seen)) << std::endl;
        return false;
    }

    void send(const std::string& keys)
    {
        if (write(master, keys.data(), keys.size()) != (ssize_t)keys.size()) throw std::runtime_error("cannot write to the pseudo terminal");
        read_for(0.1);
    }

    /**
     * Wait for the program to end and read all its output; its exit status, or -1 (also when
     * it did not end within <seconds>).
     */
    int finish(double seconds = 20)
    {
        const auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
        while (std::chrono::steady_clock::now() < end) {
            read_for(0.05);
            pollfd fd{status_fd, POLLIN, 0};
            if (poll(&fd, 1, 0) <= 0) continue;
            int code = -1;
            if (read(status_fd, &code, sizeof(code)) != (ssize_t)sizeof(code)) code = -1;
            // the program has ended, so all its output is in the terminal, which the keeper
            // still holds: read until nothing is left
            for (pollfd out{master, POLLIN, 0}; poll(&out, 1, 100) > 0;) {
                char data[4096];
                const ssize_t n = read(master, data, sizeof(data));
                if (n <= 0) break;
                output.append(data, (size_t)n);
            }
            const char byte = 0;
            if (write(release_fd, &byte, 1) != 1) kill(pid, SIGKILL);
            waitpid(pid, nullptr, 0);
            exited = true;
            return code;
        }
        return -1;
    }

    /** Whether the terminal settings are as before the program started. */
    bool restored()
    {
        termios after{};
        // the system may revoke the old descriptor when the program ends
        const int fd = open(slave_name.c_str(), O_RDWR | O_NOCTTY);
        const int r = tcgetattr(fd, &after);
        close(fd);
        const tcflag_t flags = ICANON | ECHO | ISIG;
        return r == 0 && (before.c_lflag & flags) == (after.c_lflag & flags);
    }

    std::string output;

private:
    int master = -1;
    int slave = -1;
    std::string slave_name;
    termios before{};
    int status_fd = -1;
    int release_fd = -1;
    pid_t pid = -1; // the keeper
    bool exited = false;
    size_t seen = 0;
};

/** A click at column <x>, row <y>: the SGR press and release. */
static std::string
click(int x, int y)
{
    const std::string where = std::to_string(x) + ";" + std::to_string(y);
    return "\033[<0;" + where + "M\033[<0;" + where + "m";
}

static void
test_review(const std::string& bench, const std::string& games, const std::string& work)
{
    Session session(bench, {"--live", "--review", "--tl", "--zlk", "--quick", "--sort", "--exports=extracted", games + "/vb001", games + "/vb002", games + "/vb003"}, work);
    check("review screen", session.expect("Finished."));
    check("raw input while the screen is shown", !session.restored());
    session.send("g");
    check("g sorts by name", session.expect("GAME ^"));
    // on 100 columns the game names take 25, and the first solver starts at 28
    session.send(click(30, 4));
    check("a click on a header sorts by that solver", session.expect("tl            ^"));
    session.send(click(30, 4));
    check("a second click reverses", session.expect("tl            v"));
    session.send(click(5, 6));
    session.send("e");
    check("e extracts the selected game", session.expect("Saved extracted/"));
    session.send("q");
    const int status = session.finish();
    check("q ends the review normally", status == 0);
    check("summary after the screen", session.expect("Benchmark summary"));
    const std::string& out = session.output;
    check("alternate screen entered and left", out.find("\033[?1049h") != std::string::npos && out.find("\033[?1049l") != std::string::npos);
    check("mouse reporting switched on and off", out.find("\033[?1000h\033[?1006h") != std::string::npos && out.find("\033[?1000l\033[?1006l") != std::string::npos);
    check("cursor shown again", out.rfind("\033[?25h") != std::string::npos);
    check("terminal settings restored", session.restored());
    bool extracted = false;
    for (const auto& entry : fs::recursive_directory_iterator(fs::path(work) / "extracted")) {
        if (entry.path().extension() == ".pg") extracted = true;
    }
    check("extracted game saved", extracted);
}

static void
test_interrupt(const std::string& bench, const std::string& games, const std::string& work)
{
    // a solver that takes long enough to interrupt
    Session session(bench, {"--live", "--tl", "-e", "slow:sleep 3", "--quick", "--timeout=10", games + "/vb001", games + "/vb002"}, work);
    check("live screen", session.expect("SOLVER BENCHMARK"));
    session.read_for(0.5);
    session.send("\003"); // Ctrl-C
    const int status = session.finish(8);
    check("Ctrl-C interrupts", status == 130);
    check("interrupted summary", session.expect("Interrupted"));
    check("terminal restored after Ctrl-C", session.restored());
}

static void
test_setup(const std::string& bench, const std::string& games, const std::string& work)
{
    Session session(bench, {"-i", "--tl", "--quick", "--sort", "--jsonl=setup.jsonl", games + "/vb001", games + "/vb002"}, work);
    check("setup screen", session.expect("BENCHMARK SETUP"));
    session.send("r");
    check("repetitions editor", session.expect("MEASURED REPETITIONS"));
    session.send("\0252\r"); // Ctrl-U clears the old value
    check("repetitions set", session.expect("Reps 2"));
    session.send("/zielonka recursion\r");
    check("filtered", session.expect("AVAILABLE / zielonka recursion"));
    session.send("\r");
    check("solver added", session.expect("SELECTED (2)"));
    session.send("s");
    check("benchmark runs and ends in the review", session.expect("Finished."));
    session.send("q");
    check("normal end", session.finish() == 0);
    check("summary of both solvers", session.expect("zlk"));
    check("terminal restored after setup and review", session.restored());

    // q in the setup screen ends without a benchmark
    Session cancel(bench, {"-i", "--tl", games + "/vb001"}, work);
    check("setup screen again", cancel.expect("BENCHMARK SETUP"));
    cancel.send("q");
    check("cancel", cancel.finish() == 0);
    check("no benchmark after cancel", cancel.output.find("Benchmark summary") == std::string::npos);
    check("terminal restored after cancel", cancel.restored());
}

static void
test_experiment(const std::string& bench, const std::string& oink, const std::string& games, const std::string& work)
{
    // an external solver that is slow, but correct
    const std::string slow = "slow:sleep 0.3; " + oink + " %I %O";
    std::vector<std::string> args = {"--experiment=exp", "--tl", "-e", slow, "--quick", "--sort"};
    for (const char* game : {"/vb001", "/vb002", "/vb003", "/vb004", "/vb006"}) args.push_back(games + game);
    {
        Session session(bench, args, work);
        check("experiment starts", session.expect("Experiment: "));
        session.read_for(0.8);
        session.send("\003");
        check("interrupted experiment", session.finish() == 130);
        // after the summary (the experiment also shows how to continue when it starts)
        const std::string screen = ui::strip(session.output);
        const size_t summary = screen.rfind("Interrupted; the completed results are kept.");
        check("how to continue", summary != std::string::npos && screen.find("--resume=", summary) != std::string::npos);
        check("terminal restored after the experiment", session.restored());
    }
    // --experiment again offers to resume, and ignores the other options
    Session again(bench, {"--experiment=exp", "--zlk"}, work);
    check("resume prompt", again.expect("Resume the experiment? [Y/n]"));
    again.send("\r");
    check("resumed experiment finishes", again.finish() == 0);
    check("all games finished", again.expect("(5/5 games)"));
    // the output shows paths of the test, which may contain any text, such as "zlk"
    std::string output = again.output;
    for (const std::string& path : {bench, oink, games, work}) {
        for (size_t at; (at = output.find(path)) != std::string::npos;) output.erase(at, path.size());
    }
    check("the saved solvers", output.find("zlk") == std::string::npos);
}

static void
test_timeout_review(const std::string& bench, const std::string& games, const std::string& work)
{
    Session session(bench, {"--live", "--review", "-e", "slow:sleep 5", "--repeat=3", "--warmup=0",
                            "--min-solve-time=0", "--timeout=0.3", games + "/vb001"}, work);
    check("timeout review finishes", session.expect("Finished."));
    const auto screen = ui::strip(session.output);
    check("review marks the skipped timeout", screen.find("timeout+skip") != std::string::npos &&
          screen.find("skipped timeouts") != std::string::npos);
    session.send("q");
    check("a skipped timeout still fails coverage", session.finish() == 1);
}

static void
test_timeout_resume(const std::string& bench, const std::string& games, const std::string& work)
{
    const auto root = fs::path(work) / "timeout-resume";
    const auto journal = root / "games/0.rounds";
    {
        Session session(bench, {"--fresh-process", "--experiment=timeout-resume", "-e", "slow:sleep 5",
                                "--repeat=4", "--warmup=0", "--min-solve-time=0", "--timeout=2",
                                "--skip-after-timeouts=2",
                                games + "/vb001"}, work);
        check("timeout experiment starts", session.expect("Experiment: "));
        const auto until = std::chrono::steady_clock::now() + std::chrono::seconds(15);
        bool recorded = false;
        while (std::chrono::steady_clock::now() < until) {
            std::error_code error;
            const auto size = fs::file_size(journal, error);
            if (!error && size != 0) {
                recorded = true;
                break;
            }
            session.read_for(0.05);
        }
        check("the first round finishes before interruption", recorded);
        session.send("\003");
        check("the timeout experiment is interrupted", session.finish() == 130);
    }
    // Simulate a truncated last write as well as an interrupted trial.
    std::ofstream(journal, std::ios::app) << "1 1024 1";
    Session resumed(bench, {"--resume=timeout-resume"}, work);
    check("the resumed experiment finishes with timeouts", resumed.finish() == 1);
    std::ifstream events(root / "results.jsonl");
    std::vector<std::string> trials;
    for (std::string line; std::getline(events, line);) {
        if (line.find("\"event\":\"sample\"") != std::string::npos) trials.push_back(line);
    }
    check("resume kept the completed round and counted its timeout", trials.size() == 4 &&
          trials[0].find("\"executed\":true") != std::string::npos &&
          trials[0].find("\"experiment_session\":1") != std::string::npos &&
          trials[1].find("\"executed\":true") != std::string::npos &&
          trials[1].find("\"experiment_session\":2") != std::string::npos &&
          trials[2].find("\"skip_reason\":\"timeout_limit\"") != std::string::npos &&
          trials[3].find("\"skip_reason\":\"timeout_limit\"") != std::string::npos);
    check("the journal was committed", !fs::exists(journal));
}

int
main(int argc, char** argv)
{
    if (argc != 5) {
        std::cerr << "usage: " << argv[0] << " <bench_solvers> <oink> <games> <work directory>" << std::endl;
        return 2;
    }
    signal(SIGPIPE, SIG_IGN);
    const std::string bench = argv[1];
    const std::string oink = argv[2];
    const std::string games = argv[3];
    const std::string work = argv[4];
    fs::remove_all(work);
    fs::create_directories(work);
    try {
        test_review(bench, games, work);
        test_interrupt(bench, games, work);
        test_setup(bench, games, work);
        test_experiment(bench, oink, games, work);
        test_timeout_review(bench, games, work);
        test_timeout_resume(bench, games, work);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << std::endl;
        return 1;
    }
    if (failures) std::cerr << failures << " terminal session test(s) failed" << std::endl;
    return failures ? 1 : 0;
}
