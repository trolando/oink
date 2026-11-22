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

#ifndef ZLKWJ_HPP
#define ZLKWJ_HPP

#include <cstdint>
#include <deque>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "oink/solver.hpp"

namespace pg {

/**
 * The Warsaw quasi-polynomial version of Zielonka's algorithm (as ZLKW) with the justifications
 * of Lapauw, Bruynooghe and Denecker, "Improving Parity Game Solvers with Justifications" (VMCAI
 * 2020), implicit as in ZLKJ. See zlkwj.cpp.
 */
class ZLKWJSolver : public Solver
{
public:
    ZLKWJSolver(Oink& oink, Game& game);
    virtual ~ZLKWJSolver();

    virtual void run();

protected:
    uint64_t iterations = 0;  // Zielonka steps that solve a subgame or attract the heads
    uint64_t calls = 0;       // recursive calls on a nonempty subgame
    uint64_t presumed = 0;    // steps that ran out of precision
    uint64_t reused = 0;      // calls that took the result of an earlier call
    uint64_t resets = 0;      // vertices reset because their justification was lost
    uint64_t reattracted = 0; // reset vertices attracted back

    int min_dominion = 2;     // the size of the smallest possible dominion (1 if there is a self-loop)
    int next_r = 0;           // the next region id

    /**
     * Per vertex: DIS (disabled), BOT (open), a negative escape counter (open), or the id of the
     * region that holds it. Region ids grow with time: a call takes a new id when it starts, and
     * its subgame is then the set of open vertices; afterwards it is the set of vertices with a
     * region id at least that id that are open or won by the player of the call.
     */
    std::unique_ptr<int[]> region;
    std::unique_ptr<int8_t[]> winning;
    /**
     * The strategy, which is also the justification: the successor of a vertex won by its owner,
     * -1 for a vertex won by the opponent of its owner (all its edges), and HEAD for a head of a
     * call that has no strategy yet.
     */
    std::unique_ptr<int[]> strategy;

    // A superset of the open vertices; the head search skips the other vertices with it.
    bitset open;
    // The vertices whose winner is presumed: the precision ran out where they were solved.
    bitset unjust;
    // Marks for removing duplicates from the vertices won by the player of a call.
    bitset seen;

    /**
     * The result of the last completed call: its subgame split into the vertices won by
     * Even (W[0]) and by Odd (W[1]).
     */
    std::vector<int> W[2];
    // Scratch: the vertices that are reset in a step.
    std::vector<int> reset;

    /**
     * The vectors of a call, one set per recursion depth, reused by all calls at that depth.
     * H: the heads of the call.
     * G: the vertices won by the player of the call (with stale entries of vertices that were
     *    reset or attracted by the opponent, and duplicates).
     * Y: the vertices the opponent won (the removed opponent attractors).
     * U: the pending vertices: won by the player of the call, but their justification edge was
     *    lost, or their winner is presumed. The next step that does work resets them.
     */
    struct Level {
        std::vector<int> H, G, Y, U;
    };
    std::deque<Level> levels;
    int depth = 0;

    /**
     * Results of calls in which a precision ran out, with the filter, the slots and the budget
     * of ZLKW: the vertices won by Even, then those won by Odd (a presumed vertex v as ~v), then
     * the strategies of both lists, in <V>.
     */
    static const int SLOTS = 4;
    static const size_t BUDGET = size_t(1) << 22;
    struct Prec {
        long pc, po;
        bool rbc, rbo;
    };
    struct Slot {
        uint64_t key = 0, h0 = 0;
        uint32_t n0 = 0, n1 = 0;
        std::vector<int> V;
        std::vector<Prec> precs;
        uint32_t prev = UINT32_MAX, next = UINT32_MAX; // in the list of its depth
        int depth = -1;
    };
    std::vector<Slot> slots;
    std::vector<uint32_t> depth_head, depth_tail, depth_count;
    std::vector<uint32_t> free_slots;

    /**
     * The slots of each depth form a list, least recently used first.
     */
    void unlink(uint32_t t)
    {
        Slot &M = slots[t];
        if (M.prev != UINT32_MAX) slots[M.prev].next = M.next;
        else depth_head[M.depth] = M.next;
        if (M.next != UINT32_MAX) slots[M.next].prev = M.prev;
        else depth_tail[M.depth] = M.prev;
        M.prev = M.next = UINT32_MAX;
        depth_count[M.depth]--;
    }

    void link(uint32_t t, int d)
    {
        Slot &M = slots[t];
        M.depth = d;
        M.prev = depth_tail[d];
        M.next = UINT32_MAX;
        if (M.prev != UINT32_MAX) slots[M.prev].next = t;
        else depth_head[d] = t;
        depth_tail[d] = t;
        depth_count[d]++;
    }

    void touch(uint32_t t)
    {
        const int d = slots[t].depth;
        unlink(t);
        link(t, d);
    }
    std::unordered_multimap<uint64_t, uint32_t> index;
    bool memo_on = false; // whether a precision ran out yet; only then subgames are hashed
    uint64_t res_hash[2] = {0, 0}; // after a call with a hash: the hashes of W[0] and W[1]

    /**
     * Fingerprints of calls in which a precision ran out (their subgame and precisions), in a
     * table of FILTER entries indexed by the fingerprint. A call stores its result only if its
     * fingerprint is already there, so only results of calls that came twice are stored.
     */
    static const size_t FILTER = size_t(1) << 16;
    std::vector<uint64_t> filter;

    bool seen_before(uint64_t key, long pc, long po)
    {
        const uint64_t fp = vhash(key ^ vhash((uint64_t)pc * 0x100000001b3ULL + (uint64_t)po)) | 1;
        uint64_t &e = filter[fp >> 48];
        if (e == fp) return true;
        e = fp;
        return false;
    }
    size_t pool = 0;

    static inline uint64_t vhash(uint64_t v)
    {
        // the finaliser of splitmix64
        v += 0x9e3779b97f4a7c15ULL;
        v = (v ^ (v >> 30)) * 0xbf58476d1ce4e5b9ULL;
        v = (v ^ (v >> 27)) * 0x94d049bb133111ebULL;
        return v ^ (v >> 31);
    }

    int new_region()
    {
        if (next_r == INT32_MAX) throw std::runtime_error("zlkwj: out of region ids");
        return next_r++;
    }

    int attractHeads(int i, int pl, int r, std::vector<int> &A, std::vector<int> &H);
    long attractLosing(int r, int pl, std::vector<int> &O, const std::vector<int> &H, std::vector<int> &U);
    void resetClosure(int r, int pl, std::vector<int> &U);
    int reattract(int r, int r2, int pl, std::vector<int> &A);
    bool reuse(long n, long pc, long po, uint64_t key, bool &rb_cur, bool &rb_opo);
    void store(int d, long pc, long po, uint64_t key, bool rbc, bool rbo);
    uint64_t hash_open(int top, long n);
    void solve(int top, long n, long pc, long po, bool &rb_cur, bool &rb_opo, uint64_t hg, bool hashed);
};

}

#endif
