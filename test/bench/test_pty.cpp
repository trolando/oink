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
 * Usage: test_pty <bench_solvers> <directory of test games> <work directory>
 */

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <filesystem>
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

        pid = fork();
        if (pid == 0) {
            setsid();
            const int tty = open(slave_name.c_str(), O_RDWR);
            ioctl(tty, TIOCSCTTY, 0);
            dup2(tty, 0);
            dup2(tty, 1);
            dup2(tty, 2);
            if (chdir(directory.c_str()) != 0) _exit(126);
            setenv("TERM", "xterm-256color", 1);
            unsetenv("NO_COLOR");
            std::vector<char*> argv;
            argv.push_back(const_cast<char*>(program.c_str()));
            for (const auto& arg : args) argv.push_back(const_cast<char*>(arg.c_str()));
            argv.push_back(nullptr);
            execv(program.c_str(), argv.data());
            _exit(127);
        }
    }

    ~Session()
    {
        if (pid > 0 && !exited) {
            kill(pid, SIGKILL);
            waitpid(pid, nullptr, 0);
        }
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
            if (n <= 0) return; // the program closed the terminal
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

    /** Wait for the program to end; its exit status, or -1. */
    int finish(double seconds = 20)
    {
        const auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
        while (std::chrono::steady_clock::now() < end) {
            read_for(0.05);
            int status = 0;
            if (waitpid(pid, &status, WNOHANG) == pid) {
                exited = true;
                read_for(0.1);
                return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
            }
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
    pid_t pid = -1;
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

int
main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "usage: " << argv[0] << " <bench_solvers> <games> <work directory>" << std::endl;
        return 2;
    }
    signal(SIGPIPE, SIG_IGN);
    const std::string work = argv[3];
    fs::remove_all(work);
    fs::create_directories(work);
    try {
        test_review(argv[1], argv[2], work);
        test_interrupt(argv[1], argv[2], work);
        test_setup(argv[1], argv[2], work);
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << std::endl;
        return 1;
    }
    if (failures) std::cerr << failures << " terminal session test(s) failed" << std::endl;
    return failures ? 1 : 0;
}
