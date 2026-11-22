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

#ifndef ZLKW_HPP
#define ZLKW_HPP

#include <cstdint>
#include <deque>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include "oink/solver.hpp"

namespace pg {

/**
 * The Warsaw quasi-polynomial version of Zielonka's algorithm (Parys 2019, with the
 * refinements of Lehtinen, Parys, Schewe and Wojtczak 2021, as in their implementation
 * ZLKPP-WAW), on the data structures of ZLK. See zlkw.cpp.
 */
class ZLKWSolver : public Solver
{
public:
    ZLKWSolver(Oink& oink, Game& game);
    virtual ~ZLKWSolver();

    virtual void run();

protected:
    uint64_t iterations = 0; // Zielonka steps (attractor and recursive call), as counted by ZLKPP
    uint64_t calls = 0;      // recursive calls on a nonempty subgame
    uint64_t presumed = 0;   // steps that ran out of precision
    uint64_t reused = 0;     // calls that took the result of an earlier call (see reuse())

    int min_dominion = 2;    // the size of the smallest possible dominion (1 if there is a self-loop)
    int next_r = 0;          // the next region id

    /**
     * Per vertex: DIS (disabled), BOT (open), a negative escape counter (open), or the id of the
     * region that holds it. Region ids grow with time, so a call's subgame is the set of open
     * vertices when it starts, and afterwards the vertices with a region id >= its first id.
     */
    std::unique_ptr<int[]> region;
    std::unique_ptr<int8_t[]> winning;
    std::unique_ptr<int[]> strategy;

    /**
     * The open vertices (BOT or an escape counter): those of the current subgame that are not
     * yet attracted. Exact: every attraction removes the vertex, every reset adds it.
     */
    bitset open;

    /**
     * The result of the last completed call: its subgame split into the vertices won by
     * Even (W[0]) and by Odd (W[1]).
     */
    std::vector<int> W[2];

    /**
     * The vectors of a call, one set per recursion depth, reused by all calls at that depth.
     * A: the region of the current step (the attractor of the heads), in attraction order.
     * H: the heads of A.
     * G: the subgame of the call after the last step (A and the part of the lower subgame won
     *    by the player of the call, minus what the opponent attracted).
     * Y: the vertices the opponent won so far (the removed opponent attractors).
     */
    struct Level {
        std::vector<int> A, H, G, Y;
    };
    std::deque<Level> levels;
    int depth = 0;

    /**
     * Results of calls in which a precision ran out, for a later call on the same subgame (see
     * reuse() and store() in zlkw.cpp). A slot holds the hash of a subgame and the result of a
     * call on it: the vertices won by Even, then those won by Odd, then the strategies of both
     * lists, in <V>. Several calls may have computed the same result; <precs> holds the
     * precisions of each and its flags. Every depth may have SLOTS slots, and more while all
     * slots together hold at most BUDGET ints (see store()). <free_slots> lists the unused
     * slots, <index> the slots of each hash, and <pool> counts the ints in all slots.
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

    /**
     * The hash of a vertex; the hash of a subgame is the sum over its vertices.
     */
    static inline uint64_t vhash(uint64_t v)
    {
        // the finaliser of splitmix64
        v += 0x9e3779b97f4a7c15ULL;
        v = (v ^ (v >> 30)) * 0xbf58476d1ce4e5b9ULL;
        v = (v ^ (v >> 27)) * 0x94d049bb133111ebULL;
        return v ^ (v >> 31);
    }

    bool reuse(long n, long pc, long po, uint64_t key, bool &rb_cur, bool &rb_opo);
    void store(int d, long pc, long po, uint64_t key, bool rbc, bool rbo);

    int new_region()
    {
        if (next_r == INT32_MAX) throw std::runtime_error("zlkw: out of region ids");
        return next_r++;
    }

    int attract(int i, int pl, int r, std::vector<int> &A, std::vector<int> &H);
    long attractLosing(int r, int pl, const std::vector<int> &H, std::vector<int> &R);
    uint64_t hash_open(int top, long n);
    void solve(int i, int pl, long n, long pc, long po, bool &rb_cur, bool &rb_opo, uint64_t hg, bool hashed);
};

}

#endif
