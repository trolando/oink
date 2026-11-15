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

#include "worker.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <thread>

#ifndef _WIN32
#include <cerrno>
#include <fcntl.h>
#include <poll.h>
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace bench {

#ifndef _WIN32

namespace {

using Clock = std::chrono::steady_clock;

// The end of the child's output kept for the message of a crash.
const size_t MAX_OUTPUT = 16384;

/** A request from the benchmark to the child: run one trial. */
struct Request
{
    uint64_t solver;
    uint32_t seed;
};

/** The fixed part of a reply; the three strings of the trial follow it. */
struct ReplyHeader
{
    int32_t status;
    int32_t parallel;
    int32_t workers;
    uint32_t sizes[3];
    double times[4];
};

bool
write_all(int fd, const void* data, size_t size)
{
    const char* p = (const char*)data;
    while (size > 0) {
        const ssize_t n = ::write(fd, p, size);
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) return false;
        p += n;
        size -= (size_t)n;
    }
    return true;
}

bool
read_all(int fd, void* data, size_t size)
{
    char* p = (char*)data;
    while (size > 0) {
        const ssize_t n = ::read(fd, p, size);
        if (n < 0 && errno == EINTR) continue;
        if (n <= 0) return false;
        p += n;
        size -= (size_t)n;
    }
    return true;
}

void
send_reply(int fd, const Trial& trial)
{
    ReplyHeader header;
    std::memset(&header, 0, sizeof(header));
    header.status = (int32_t)trial.status;
    header.parallel = trial.parallel ? 1 : 0;
    header.workers = trial.workers;
    header.sizes[0] = (uint32_t)trial.effective_solver.size();
    header.sizes[1] = (uint32_t)trial.effective_config.size();
    header.sizes[2] = (uint32_t)trial.message.size();
    header.times[0] = trial.setup;
    header.times[1] = trial.solve;
    header.times[2] = trial.verify;
    header.times[3] = trial.cpu;
    std::string data((const char*)&header, sizeof(header));
    data += trial.effective_solver;
    data += trial.effective_config;
    data += trial.message;
    if (!write_all(fd, data.data(), data.size())) _exit(4);
}

/** Decode a complete reply from <data>; false if more bytes are needed. */
bool
parse_reply(const std::string& data, Trial& trial)
{
    if (data.size() < sizeof(ReplyHeader)) return false;
    ReplyHeader header;
    std::memcpy(&header, data.data(), sizeof(header));
    const size_t total = sizeof(header) + (size_t)header.sizes[0] + header.sizes[1] + header.sizes[2];
    if (data.size() < total) return false;
    trial.status = (Status)header.status;
    trial.parallel = header.parallel != 0;
    trial.workers = header.workers;
    trial.setup = header.times[0];
    trial.solve = header.times[1];
    trial.verify = header.times[2];
    trial.cpu = header.times[3];
    size_t at = sizeof(header);
    trial.effective_solver = data.substr(at, header.sizes[0]);
    at += header.sizes[0];
    trial.effective_config = data.substr(at, header.sizes[1]);
    at += header.sizes[1];
    trial.message = data.substr(at, header.sizes[2]);
    return true;
}

void
close_fd(int& fd)
{
    if (fd != -1) ::close(fd);
    fd = -1;
}

}

bool
Worker::supported()
{
    return true;
}

Worker::Worker(Body body, Hook begin, Hook end) : body(std::move(body)), begin(std::move(begin)), end(std::move(end))
{
}

Worker::~Worker()
{
    stop();
}

bool
Worker::running() const
{
    return pid > 0;
}

void
Worker::start()
{
    int requests[2];
    int replies[2];
    int out[2];
    if (pipe(requests) != 0) throw std::runtime_error("cannot create a pipe for the worker");
    if (pipe(replies) != 0) {
        ::close(requests[0]);
        ::close(requests[1]);
        throw std::runtime_error("cannot create a pipe for the worker");
    }
    if (pipe(out) != 0) {
        for (int fd : {requests[0], requests[1], replies[0], replies[1]}) ::close(fd);
        throw std::runtime_error("cannot create a pipe for the worker");
    }
    // external solvers started by the child must not hold these pipes
    for (int fd : {requests[0], requests[1], replies[0], replies[1], out[0], out[1]}) fcntl(fd, F_SETFD, FD_CLOEXEC);

    // the child gets a copy of everything buffered: flush it first
    std::cout.flush();
    std::cerr.flush();
    std::fflush(nullptr);

    const pid_t child = fork();
    if (child < 0) {
        for (int fd : {requests[0], requests[1], replies[0], replies[1], out[0], out[1]}) ::close(fd);
        throw std::runtime_error("cannot start a worker process");
    }

    if (child == 0) {
        // The worker: its own process group, default signals, no terminal input.
        setpgid(0, 0);
        signal(SIGINT, SIG_DFL);
        signal(SIGTERM, SIG_DFL);
        signal(SIGPIPE, SIG_DFL);
        const int null = open("/dev/null", O_RDONLY);
        if (null >= 0) dup2(null, STDIN_FILENO);
        dup2(out[1], STDOUT_FILENO);
        dup2(out[1], STDERR_FILENO);
        ::close(requests[1]);
        ::close(replies[0]);
        ::close(out[0]);
        try {
            if (begin) begin();
            Request request;
            while (read_all(requests[0], &request, sizeof(request))) {
                const Trial trial = body((size_t)request.solver, request.seed);
                send_reply(replies[1], trial);
            }
            if (end) end();
        } catch (const std::exception& error) {
            std::fprintf(stderr, "worker error: %s\n", error.what());
            _exit(3);
        }
        // no destructors or exit handlers: they belong to the benchmark process
        _exit(0);
    }

    setpgid(child, child); // also here, in case the child has not done it yet
    ::close(requests[0]);
    ::close(replies[1]);
    ::close(out[1]);
    pid = child;
    request_fd = requests[1];
    reply_fd = replies[0];
    output_fd = out[0];
    fcntl(reply_fd, F_SETFL, O_NONBLOCK);
    fcntl(output_fd, F_SETFL, O_NONBLOCK);
    reply.clear();
    output.clear();
    generation_++;
}

void
Worker::drain_output()
{
    char data[4096];
    for (;;) {
        const ssize_t n = ::read(output_fd, data, sizeof(data));
        if (n <= 0) break;
        output.append(data, (size_t)n);
    }
    if (output.size() > MAX_OUTPUT) output.erase(0, output.size() - MAX_OUTPUT);
}

/**
 * End the child: <force> kills its process group at once; otherwise the child
 * gets a moment to finish after its requests are closed.
 */
void
Worker::kill(bool force)
{
    if (pid <= 0) return;
    close_fd(request_fd);
    bool reaped = false;
    if (!force) {
        for (int i = 0; i < 100 && !reaped; i++) {
            int status = 0;
            const pid_t r = waitpid((pid_t)pid, &status, WNOHANG);
            if (r == (pid_t)pid || (r < 0 && errno != EINTR)) {
                reaped = true;
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(2));
            }
        }
    }
    if (!reaped) {
        // the group includes the processes started by external solvers
        ::kill(-(pid_t)pid, SIGKILL);
        ::kill((pid_t)pid, SIGKILL);
        int status = 0;
        while (waitpid((pid_t)pid, &status, 0) < 0 && errno == EINTR) {
        }
    }
    close_fd(reply_fd);
    close_fd(output_fd);
    pid = -1;
}

void
Worker::stop()
{
    kill(false);
}

/** The child ended without a reply: collect its status and output. */
Trial
Worker::crashed()
{
    Trial trial;
    trial.status = Status::Crash;
    int status = 0;
    pid_t r = 0;
    for (int i = 0; i < 500; i++) {
        r = waitpid((pid_t)pid, &status, WNOHANG);
        if (r != 0) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    drain_output();
    if (r == (pid_t)pid && WIFSIGNALED(status)) {
        const int sig = WTERMSIG(status);
        trial.message = "worker terminated by signal " + std::to_string(sig) + " (" + strsignal(sig) + ")";
    } else if (r == (pid_t)pid && WIFEXITED(status)) {
        trial.message = "worker exited with status " + std::to_string(WEXITSTATUS(status)) + " without a result";
    } else {
        trial.message = "worker closed its output without a result";
    }
    if (!output.empty()) trial.message += "\n" + output;
    if (r == (pid_t)pid) {
        // reaped: only close our ends, the process group is gone with it
        close_fd(request_fd);
        close_fd(reply_fd);
        close_fd(output_fd);
        pid = -1;
    } else {
        kill(true);
    }
    return trial;
}

Trial
Worker::run(size_t solver, uint32_t seed, double timeout, volatile std::sig_atomic_t* quit)
{
    if (!running()) start();
    reply.clear();
    Request request;
    std::memset(&request, 0, sizeof(request));
    request.solver = solver;
    request.seed = seed;
    if (!write_all(request_fd, &request, sizeof(request))) return crashed();

    const auto begin_time = Clock::now();
    for (;;) {
        if (quit != nullptr && *quit) {
            kill(true);
            Trial trial;
            trial.status = Status::Interrupted;
            trial.message = "interrupted";
            return trial;
        }
        const double elapsed = std::chrono::duration<double>(Clock::now() - begin_time).count();
        if (timeout > 0 && elapsed >= timeout) {
            kill(true);
            Trial trial;
            trial.status = Status::Timeout;
            trial.message = "the trial exceeded the wall-clock limit (setup, solve and verification)";
            return trial;
        }
        // wake up regularly to check the limit and the interrupt flag
        double wait = 0.2;
        if (timeout > 0 && timeout - elapsed < wait) wait = timeout - elapsed;
        pollfd fds[2] = {{reply_fd, POLLIN, 0}, {output_fd, POLLIN, 0}};
        const int ready = poll(fds, 2, (int)(wait * 1000) + 1);
        if (ready < 0 && errno != EINTR) return crashed();
        if (ready <= 0) continue;
        if (fds[1].revents) drain_output();
        if (fds[0].revents) {
            char data[65536];
            bool closed = false;
            for (;;) {
                const ssize_t n = ::read(reply_fd, data, sizeof(data));
                if (n > 0) {
                    reply.append(data, (size_t)n);
                    continue;
                }
                if (n == 0) closed = true;
                break;
            }
            Trial trial;
            if (parse_reply(reply, trial)) {
                drain_output();
                output.clear();
                return trial;
            }
            if (closed) return crashed();
        }
    }
}

#else

bool
Worker::supported()
{
    return false;
}

Worker::Worker(Body body, Hook begin, Hook end) : body(std::move(body)), begin(std::move(begin)), end(std::move(end))
{
}

Worker::~Worker()
{
}

bool
Worker::running() const
{
    return false;
}

void
Worker::start()
{
    throw std::runtime_error("worker processes are not available on this platform");
}

void
Worker::kill(bool)
{
}

void
Worker::stop()
{
}

void
Worker::drain_output()
{
}

Trial
Worker::crashed()
{
    return Trial();
}

Trial
Worker::run(size_t, uint32_t, double, volatile std::sig_atomic_t*)
{
    Trial trial;
    trial.status = Status::SetupError;
    trial.message = "worker processes are not available on this platform; use --in-process";
    return trial;
}

#endif

}
