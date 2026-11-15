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
 * Process isolation of the trials.
 *
 * A worker is a child process that runs the trials of one game, so a crash or
 * a hang of a solver cannot take the benchmark down, and a timeout can stop
 * it. The child is forked from the benchmark, so it already has the game in
 * memory. It serves the trials of its game in turn, sharing the warm caches
 * and allocator state between them as one process would, until it is stopped.
 * After a failure the benchmark starts a fresh worker.
 *
 * The child runs in its own process group, which also holds the processes of
 * external solvers, so a timeout stops all of them, and SIGINT from the
 * terminal reaches only the benchmark.
 *
 * Only on POSIX systems; elsewhere supported() is false and the trials run in
 * the benchmark process itself, without timeouts.
 */

#ifndef OINK_BENCH_WORKER_HPP
#define OINK_BENCH_WORKER_HPP

#include <csignal>
#include <cstdint>
#include <functional>
#include <string>

#include "bench.hpp"

namespace bench {

class Worker
{
public:
    using Body = std::function<Trial(size_t solver, uint32_t seed)>;
    using Hook = std::function<void()>;

    /** Whether this platform can isolate trials in a worker process. */
    static bool supported();

    /**
     * A worker that runs <body> for every trial. <begin> and <end> run in
     * every child, before its first and after its last trial (to start and
     * stop Lace, for example); either may be empty.
     */
    Worker(Body body, Hook begin = Hook(), Hook end = Hook());
    ~Worker();
    Worker(const Worker&) = delete;
    Worker& operator=(const Worker&) = delete;

    /**
     * Run the trial of solver <solver> with seed <seed> in the child, starting
     * a child if none runs. Waits at most <timeout> seconds (0: no limit) and
     * stops early when <*quit> is set (if given).
     */
    Trial run(size_t solver, uint32_t seed, double timeout, volatile std::sig_atomic_t* quit);

    /** Stop the child and every process in its group. */
    void stop();

    bool running() const;

    /** The number of children started so far. */
    int generation() const { return generation_; }

private:
    Body body;
    Hook begin;
    Hook end;
    int generation_ = 0;
    long pid = -1;
    int request_fd = -1;
    int reply_fd = -1;
    int output_fd = -1;
    std::string reply;      // the reply received so far
    std::string output;     // the end of what the child wrote to stdout and stderr

    void start();
    void kill(bool force);
    void drain_output();
    Trial crashed();
};

}

#endif
