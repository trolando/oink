/*
 * Copyright 2022-2026 Tom van Dijk
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
 *
 * The Strahler measures are adapted from the STRPM solver by Vanessa Flügel and
 * Guillermo A. Pérez (https://github.com/vfluegel/oink), see gpm.hpp.
 */

#include <algorithm>
#include <cassert>
#include <cstring>

#include "gpm.hpp"
#include "oink/uintqueue.hpp"

namespace pg {


// Append an integral value T to buf in little-endian.
template <class T>
inline void push_le(std::vector<std::byte>& buf, T v) {
    static_assert(std::is_integral_v<T>, "push_le: T must be an integral type");
    using U = std::make_unsigned_t<T>;
    U x;
    static_assert(sizeof(x) == sizeof(v), "sizes must match");
    // bitwise copy to avoid implementation-defined casts for signed types
    std::memcpy(&x, &v, sizeof(T));
    for (size_t i = 0; i < sizeof(T); ++i) {
        buf.push_back(std::byte((x >> (8*i)) & U{0xFF}));
    }
}

// Read an integral value T from buf[offset..] (little-endian).
template <class T>
inline T read_le_at(const std::vector<std::byte>& buf, size_t offset) {
    static_assert(std::is_integral_v<T>, "read_le_at: T must be an integral type");
    using U = std::make_unsigned_t<T>;
    U x = 0;
    for (size_t i = 0; i < sizeof(T); ++i) {
        x |= U(std::to_integer<uint8_t>(buf[offset + i])) << (8*i);
    }
    T v;
    static_assert(sizeof(v) == sizeof(x), "sizes must match");
    std::memcpy(&v, &x, sizeof(T));  // preserves bit pattern (handles signed)
    return v;
}

template <class T>
inline T read_le_advance(const std::vector<std::byte>& buf, size_t& cursor) {
    T v = read_le_at<T>(buf, cursor);
    cursor += sizeof(T);
    return v;
}


Small::Small(const Game &game, const int player) : game(game)
{
    this->player = player;
    this->n = game.nodecount();                  // n := number of vertices
    this->d = game.priority(game.nodecount()-1); // d := highest priority
    this->l = (d+2)/2;                           // l := length of the tuples

    // count how many vertices of each priority there are
    this->counts = new int[d+1];
    for (int p=0; p<=d; p++) this->counts[p] = 0;
    for (int v=0; v<n; v++) this->counts[game.priority(v)]++;

    // allocate and initialize the vertices
    data = new int[l * (n+4)];
    std::fill(data, data + (l*(n+4)), '\0');
    top.resize(n+4);
}

Small::Small(const Small &src) : game(src.game)
{
    this->player = src.player;
    this->n = src.n;
    this->d = src.d;
    this->l = src.l;
    this->counts = new int[d+1];
    this->data = new int[l * (n+4)];
    std::copy(src.counts, src.counts+(d+1), this->counts);
    std::copy(src.data, src.data+(l*(n+4)), this->data);
    this->top = src.top;
}

Small::~Small()
{
    delete[] data;
    delete[] counts;
}

void
Small::set(const Measures &other)
{
    auto src = dynamic_cast<const Small*>(&other);

    delete[] data;
    delete[] counts;

    this->player = src->player;
    this->n = src->n;
    this->d = src->d;
    this->l = src->l;
    this->counts = new int[d+1];
    this->data = new int[l * (n+4)];
    std::copy(src->counts, src->counts+(d+1), this->counts);
    std::copy(src->data, src->data+(l*(n+4)), this->data);
    this->top = src->top;
}

void Small::copy(int from, int to)
{
    // decrease counts?
    if (to >= 0 && top[to+4] == 0 && top[from+4] == 1) {
        int pr = game.priority(to);
        if ((pr&1)==player) {
            counts[pr]--;
        }
    }

    int *fp = data + (from+4)*l;
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = *fp++;
    top[to+4] = top[from+4];
}

void Small::copy_from(const Measures& mother, int from, int to)
{
    const Small *other = dynamic_cast<const Small*>(&mother);

    // decrease counts?
    if (to >= 0 && top[to+4] == 0 && other->top[from+4] == 1) {
        int pr = game.priority(to);
        if ((pr&1)==player) {
            counts[pr]--;
        }
    }

    int *fp = other->data + (from+4)*l;
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = *fp++;
    top[to+4] = other->top[from+4];
}

void Small::copy_bot(int to)
{
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = 0;
    top[to+4] = 0;
}

void Small::copy_top(int to)
{
    // decrease counts?
    if (to >= 0 && top[to+4] == 0) {
        int pr = game.priority(to);
        if ((pr&1)==player) {
            counts[pr]--;
        }
    }
            
    top[to+4] = 1;
}

bool Small::is_top(int tgt) const
{
    return top[tgt+4];
}

void Small::see(int tgt, int priority)
{
    if (top[tgt+4]) return; // don't care if already at top

    // for player 0 tuple of 0, 2, 4, 6, 8...
    // for player 1 tuple of 1, 3, 5, 7, 9...

    int *tp = data + (tgt+4)*l;
    for (int i=0; i<l; i++) {
        int tuple_p = this->player + 2*i;
        if (tuple_p < priority) {
            tp[i] = 0;
        } else if (tuple_p == priority) {
            // check if reached the max! if so, carry on
            if (tp[i] >= counts[tuple_p]) {
                tp[i] = 0;
                // on carry overflow, this is a \Top
                if ((priority+2) > d) {
                    top[tgt + 4] = 1;
                    return; // done
                } else {
                    priority+=2;
                }
            } else {
                tp[i]++;
                return; // done
            }
        } else {
            break;
        }
    }
}

void Small::inc(int tgt)
{
    if (top[tgt+4]) return; // nothing to do

    int *tp = data + (tgt+4)*l;
    for (int i=0; i<l; i++) {
        int tuple_p = this->player + 2*i;
        if (tp[i] == counts[tuple_p]) {
            tp[i] = 0;
        } else {
            tp[i]++;
            return;
        }
    }

    top[tgt+4] = true;
}

/**
 * Return -1 if value[i] < value[j]
 * Return  0 if value[i] = value[j]
 * Return  1 if value[i] > value[j]
 */
int Small::compare(int i, int j) const
{
    // Handle \Top cases
    if (top[i+4]) {
        if (top[j+4]) return 0;
        return 1;
    }
    if (top[j+4]) return -1;
    // Compare tuples, starting at the highest (last) value
    int *ip = data + (i+4)*l;
    int *jp = data + (j+4)*l;
    for (int x=l-1; x>=0; x--) {
        if (ip[x] < jp[x]) return -1;
        if (ip[x] > jp[x]) return 1;
    }
    // No difference? then equal
    return 0;
}

bool Small::eq(int i, int j) const
{
    if (top[i+4]) {
        if (top[j+4]) return true;
        return false;
    }
    if (top[j+4]) return false;


    int *ip = data + (i+4)*l;
    int *jp = data + (j+4)*l;
    for (int x = 0; x < l; x++) {
        if (ip[x] != jp[x]) return false;
    }
    return true;
}

void Small::stream(std::ostream &out, int i) const
{
    int *ip = data + (i+4)*l;
    if (top[i+4]) {
        out << "Top";
    } else {
        out << "[";
        for (int x=l-1; x >= 0; x--) out << " " << ip[x];
        out << " ]";
    }
}

std::vector<std::byte> Small::store(int tgt) const
{
    std::vector<std::byte> blob;
    blob.reserve(static_cast<size_t>(l) * sizeof(int) + 1); // perf: avoid reallocs

    int *tp = data + (tgt+4) * l;
    for (int i=0; i<l; i++) {
        push_le<int>(blob, *tp++);
    }
    push_le<uint8_t>(blob, top[tgt+4]);

    return blob;
}

void Small::load(int tgt, const std::vector<std::byte>& blob)
{
    size_t cur = 0;

    int* tp = data + (tgt + 4) * l;
    for (int i = 0; i < l; ++i) {
        *tp++ = read_le_advance<int>(blob, cur);
    }

    top[tgt + 4] = read_le_advance<uint8_t>(blob, cur);
}

Ordered::Ordered(const Game &game, const int player)
{
    this->player = player;
    this->n = game.nodecount();                  // n := number of vertices
    this->d = game.priority(game.nodecount()-1); // d := highest priority
    priorities.resize(d + 1);

    // find highest own priority
    this->max_pr = -1;
    this->max_opp_pr = -1;
    for (int i=game.nodecount()-1; i>=0; i--) {
        const int pr = game.priority(i);
        if ((pr&1) == player ) {
            if (pr > this->max_pr) {
                this->max_pr = pr;
            }
        } else if(pr > max_opp_pr) {
            this->max_opp_pr = pr;
        }
        priorities[pr] = true;
    }
    if (max_opp_pr == -1) {
        max_opp_pr = player;
    }

    // compute number of vertices with player's priority
    this->N = 0;
    for (auto v=0; v<game.nodecount(); v++) {
        if ((game.priority(v)&1) == player) N++;
    }

    // see the paper: use l = ceil(log2(N+1)) + 1
    l = 1;
    while ((1L<<(l-1)) <= N) l++;

    //std::cout << "n=" << N << std::endl;
    //std::cout << "l=" << l << std::endl;

    // count how many vertices of each priority there are
    this->counts = new int[d+1];
    for (int p=0; p<=d; p++) this->counts[p] = 0;
    for (int v=0; v<n; v++) this->counts[game.priority(v)]++;

    // allocate and initialize the vertices
    data = new int[l * (n+4)];
    std::fill(data, data + (l*(n+4)), -1);
    top.resize(n+4);
}

Ordered::Ordered(const Ordered &src)
{
    this->player = src.player;
    this->n = src.n;
    this->N = src.N;
    this->d = src.d;
    this->l = src.l;
    this->max_pr = src.max_pr;
    this->max_opp_pr = src.max_opp_pr;
    this->counts = new int[d+1];
    this->data = new int[l * (n+4)];
    std::copy(src.counts, src.counts+(d+1), this->counts);
    std::copy(src.data, src.data+(l*(n+4)), this->data);
    this->top = src.top;
    this->priorities = src.priorities;
}

Ordered::~Ordered()
{
    delete[] data;
    delete[] counts;
}

void
Ordered::set(const Measures &other)
{
    auto src = dynamic_cast<const Ordered*>(&other);

    delete[] data;
    delete[] counts;

    this->player = src->player;
    this->n = src->n;
    this->N = src->N;
    this->d = src->d;
    this->l = src->l;
    this->max_pr = src->max_pr;
    this->max_opp_pr = src->max_opp_pr;
    this->counts = new int[d+1];
    this->data = new int[l * (n+4)];
    std::copy(src->counts, src->counts+(d+1), this->counts);
    std::copy(src->data, src->data+(l*(n+4)), this->data);
    this->top = src->top;
    this->priorities = src->priorities;
}

void Ordered::copy(int from, int to)
{
    int *fp = data + (from+4)*l;
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = *fp++;
    top[to+4] = top[from+4];
}

void Ordered::copy_from(const Measures& mother, int from, int to)
{
    const Ordered *other = dynamic_cast<const Ordered*>(&mother);
    int *fp = other->data + (from+4)*l;
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = *fp++;
    top[to+4] = other->top[from+4];
}

void Ordered::copy_bot(int to)
{
    int *tp = data + (to+4)*l;
    for (int i=0; i<l; i++) *tp++ = -1;
    top[to+4] = 0;
}

void Ordered::copy_top(int to)
{
    top[to+4] = 1;
}

bool Ordered::is_top(int tgt) const
{
    return top[tgt+4];
}

int Ordered::val(int tgt) const
{
    int res = 0;
    int *tp = data + (tgt+4)*l;
    bool found = false;
    for (int x=l-1; x>=0; x--) {
        res *= 2;
        if (!found) {
            res += (tp[x] != -1) ? 1 : 0;
            if ((tp[x] & 1) == player) {
                found = true;
            }
        }
    }
    return res;
}

void Ordered::see(int tgt, int priority)
{
    if (top[tgt+4]) return; // nothing to do
    assert(tgt != -3); // (scratch slot; the solvers use -1)
    // the result is the smaller of <tgt> raised by <priority> (computed in place) and of
    // <tgt> first bumped and then raised (computed in -3; bump fails if the lowest value is _)
    if (data[(tgt+4)*l] == -1) {
        up(tgt, priority);
    } else {
        copy(tgt, -3);
        up(tgt, priority);
        if (bump(-3)) {
            up(-3, priority);
            if (compare(-3, tgt) < 0) {
                copy(-3, tgt);
            }
        }
    }
    int last = *(data + (tgt+4)*l + l - 1);
    if (last != -1 && (last & 1) == player) {
        top[tgt+4] = true;
    }
}

bool Ordered::up(int tgt, int p)
{
    int *tp = data + (tgt+4)*l;

    // return if already top
    if (top[tgt+4]) return false;

    // first set the first _ or non-<player> value to "priority"
    for (int i=0; i<l; i++) {
        if (tp[i] == -1 or ((tp[i] & 1) != player)) {
            tp[i] = p; // set it
            break;
        } else {
            tp[i] = -1;
        }

    }

    // next set the highest non-_ < p value to "p" and the rest to -1
    for (int i=l-1; i>=0; i--) {
        if (tp[i] == -1) continue;
        if ((p&1) == player) {
            if (tp[i] >= p) continue;
        } else {
            if (tp[i] > p) continue;
        }
        tp[i] = p;
        for (int j=i-1; j>=0; j--) tp[j] = -1; // set rest to -1
        return true;
    }

    return true;
}

bool Ordered::bump(int tgt)
{
    if (top[tgt+4]) return false; // nothing to do

    int *tp = data + (tgt+4)*l;

    if (tp[0] == -1) return false; // nothing to do

    for (int i = 1; i<l; i++) {
        if (tp[i] == -1) {
            // if the current is _, then the smallest increase is max_opp_pr
            tp[i] = max_opp_pr;
            for (int z=i-1; z>=0; z--) tp[z] = -1;
            return true;
        } else if ((tp[i] & 1) != player) {
            // if the current is of opponent parity, then subtract 2 OR change parity
            if (tp[i] > (1-player)) {
                // if odd != -1 / even != 0 then we can just subtract 2 and we are done
                tp[i] -= 2;
                for (int z=i-1; z>=0; z--) tp[z] = -1;
                return true;
            } else if (i < l-1 and tp[i+1] >= player) {
                // if is 1/0, then can set to 0/1 if the previous is set and we're not the highest...
                tp[i] = player;
                for (int z=i-1; z>=0; z--) tp[z] = -1;
                return true;
            }
        } else {
            // if it's of parity <player>, and the next one is higher, then
            // increase and done probably the only one that is relevant, since
            // this algorithm is used when an "even" sequence must be bumped
            if (tp[i] + 2 <= max_pr and i < l-1 and (tp[i+1] == -1 || tp[i+1] >= tp[i] + 2)) {
                tp[i] += 2;
                for (int z=i-1; z>=0; z--) tp[z] = -1;
                return true;
            }
        }
    }

    return false;
}

void Ordered::inc(int tgt)
{
    if (top[tgt+4]) return;

    int *tp = data + (tgt+4)*l;

    for (int i=0; i<l; i++) {
        int min = d;
        for (int j=i+1; j<l; j++) {
            if (tp[j] == -1) continue;
            if (tp[j] < min) min = tp[j];
        }
        if (tp[i] == -1) {
            tp[i] = (min&1) == player && min != 0 ? min-1 : min;
            return;
        } else if ((tp[i]&1) == player) {
            // only bump if the next higher allows it
            if ((tp[i]+2) <= min) {
                tp[i] += 2;
                return;
            }
        } else {
            // 5->3, 4->2, 3->1, 2->0, 1->0, 0->1
            if (tp[i] >= 2) {
                tp[i] -= 2;
                return;
            } else if (tp[i] == 1) {
                tp[i] = 0;
                return;
            } else if (min >= 1) {
                tp[i] = 1;
                return;
            }
        }
        tp[i] = -1;
    }

    top[tgt+4] = true; // cannot go higher, therefore top
}

/**
 * Return -1 if value[i] < value[j]
 * Return  0 if value[i] = value[j]
 * Return  1 if value[i] > value[j]
 */
int Ordered::compare(int i, int j) const
{
    // Handle \Top cases
    if (top[i+4]) {
        if (top[j+4]) return 0;
        return 1;
    }
    if (top[j+4]) return -1;

    // Compare tuples, starting at the highest (last) value
    int *ip = data + (i+4)*l;
    int *jp = data + (j+4)*l;
    for (int x=l-1; x>=0; x--) {
        if (ip[x] == jp[x]) continue;
        // cases where i.x is _ or j.x is _
        if (ip[x] == -1) return -1;
        if (jp[x] == -1) return 1;
        // case distinction based on parity
        int pr_i = ip[x];
        if ((pr_i&1) != player) pr_i = -pr_i;
        int pr_j = jp[x];
        if ((pr_j&1) != player) pr_j = -pr_j;
        if (pr_i < pr_j) return -1;
        if (pr_i > pr_j) return 1;
    }
    // No difference? then equal
    return 0;
}

bool Ordered::eq(int i, int j) const
{
    if (top[i+4]) {
        if (top[j+4]) return true;
        return false;
    }
    if (top[j+4]) return false;


    int *ip = data + (i+4)*l;
    int *jp = data + (j+4)*l;
    for (int x = 0; x < l; x++) {
        if (ip[x] != jp[x]) return false;
    }
    return true;
}

void Ordered::stream(std::ostream &out, int i) const
{
    int *ip = data + (i+4)*l;
    if (top[i+4]) {
        out << "Top";
    } else {
        out << "[";
        for (int x = l - 1; x >= 0; x--) {
            if (ip[x] == -1) out << " _";
            else out << " " << ip[x];
        }
        out << " ]";
    }
}

std::vector<std::byte> Ordered::store(int tgt) const
{
    std::vector<std::byte> blob;
    blob.reserve(static_cast<size_t>(l) * sizeof(int) + 1); // perf: avoid reallocs
        
    int *tp = data + (tgt+4) * l; 
    for (int i=0; i<l; i++) {
        push_le<int>(blob, *tp++);
    }
    push_le<uint8_t>(blob, top[tgt+4]);

    return blob;
}       
    
void Ordered::load(int tgt, const std::vector<std::byte>& blob)
{       
    size_t cur = 0;

    int* tp = data + (tgt + 4) * l;
    for (int i = 0; i < l; ++i) {
        *tp++ = read_le_advance<int>(blob, cur);
    }

    top[tgt + 4] = read_le_advance<uint8_t>(blob, cur);
}


Succinct::Succinct(const Game &game, const int player)
{
    this->player = player;
    this->N = game.nodecount();
    this->l = 1;
    while ((1L<<(this->l)) <= N) this->l++;
    int max_pr = game.priority(game.nodecount()-1);
    this->h = player == 0 ? (max_pr/2)+1 : (max_pr+1)/2;
    code.assign((size_t)l*(4+N), 0); // all bits 0 at depth 0
}

void
Succinct::set(const Measures &other)
{
    // (not via operator=: Measures::operator= calls set)
    auto src = dynamic_cast<const Succinct*>(&other);
    this->player = src->player;
    this->N = src->N;
    this->l = src->l;
    this->h = src->h;
    this->code = src->code;
}

void Succinct::copy(int from, int to)
{
    std::memcpy(&code[(size_t)l*(4+to)], &code[(size_t)l*(4+from)], sizeof(uint32_t)*l);
}

void Succinct::copy_from(const Measures& mother, int from, int to)
{
    const Succinct *other = dynamic_cast<const Succinct*>(&mother);
    // assuming that the two Measures have the same value for l...
    std::memcpy(&code[(size_t)l*(4+to)], &other->code[(size_t)l*(4+from)], sizeof(uint32_t)*l);
}

void Succinct::copy_bot(int to)
{
    std::fill_n(&code[(size_t)l*(4+to)], l, 0u);
}

void Succinct::copy_top(int to)
{
    code[(size_t)l*(4+to)] = TOP;
}

bool Succinct::is_top(int target) const
{
    return code[(size_t)l*(4+target)] == TOP;
}

// the code for bit <b> at depth <d>
static inline uint32_t succ_code(int d, uint32_t b) { return ((uint32_t)d << 1) | b; }

void Succinct::trunc(int target, int pindex)
{
    uint32_t *c = &code[(size_t)l*(4+target)];
    // check if Top
    if (c[0] != TOP) {
        // [pindex],....,....,.... => [pindex],0*
        for (int i=l-1; i>=0 and (int)(c[i]>>1) > pindex; i--) c[i] = succ_code(pindex+1, 0);
    }
}

void Succinct::see(int target, int priority)
{
    int pindex = player == 0 ? h-(priority+1)/2-1 : h-priority/2-1;

    if ((priority&1) != player) {
        trunc(target, pindex);
        return;
    }

    uint32_t *c = &code[(size_t)l*(4+target)];
    auto depth = [&](int i) { return (int)(c[i] >> 1); };
    auto bit = [&](int i) { return c[i] & 1u; };

    // Case 1: Top -> Top
    if (c[0] == TOP) return; // nothing to do, it's top!

    // Case 2: Some bits below [pindex], so we can let [pindex] go from ..e to ..10*
    if (depth(l-1) > pindex) {
        // we know the last one is below pindex. So pindex ends with epsilon currently.
        // i.e., xxxe,...,... => xxx10*
        int i;
        for (i=l-1; i>=0 and depth(i) > pindex; i--) c[i] = succ_code(pindex, 0);
        c[i+1] |= 1u;
        return;
    }

    // Case 3: no bits below [pindex], so analyze lowest nonempty level
    //         that is: all "depth" are <= pindex
    //         that is: we simply find the next higher measure
    //
    // Case 3a: the lowest has a 0 ==> to epsilon, but special case for bottom
    // ?0111 b  => ?e?
    // ?0       => ?,0
    // ?0111    => ?,0000
    // when we find 0 on deepest, set tail to 0 and deepest+1
    //
    // Case 3b: the bottom one has only 1s (cannot be empty)
    //          next one is not empty
    // ?,111 => ?100
    // set to 1, then only 0s, at deepest-1

    // Case 3c: the bottom one has only 1s (cannot be empty)
    //          next one is empty
    // ?,,111   => ?,100
    // set to 1, then only 0s, at deepest-1

    // first zero: ,00000

    int deepest = depth(l-1);
    for (int i=l-1; i>=0; i--) {
        if (depth(i) == deepest) {
            if (bit(i) == 0) {
                for (int j=i; j<l; j++) c[j] = succ_code(deepest+1, bit(j));
                return;
            } else {
                // set every later bit to 0
                c[i] &= ~1u;
            }
        } else {
            c[i+1] |= 1u;
            for (int j=i+1; j<l; j++) c[j] = succ_code(deepest-1, bit(j));
            return;
        }
    }

    // if we're here, then...
    // ,1111 => 1000,
    // 1111 => T
    if (deepest == 0) {
        c[0] = TOP; // to Top!
    } else {
        // everything else is set to 0 already
        c[0] = succ_code(deepest-1, 1);
        for (int i=1; i<l; i++) c[i] = succ_code(deepest-1, bit(i));
    }
}

void Succinct::inc(int target)
{
    see(target, player);
}

/**
 * Return -1 if value[i] < value[j]
 * Return  0 if value[i] = value[j]
 * Return  1 if value[i] > value[j]
 */
int Succinct::compare(int u, int v) const
{
    const uint32_t *cu = &code[(size_t)l*(4+u)];
    const uint32_t *cv = &code[(size_t)l*(4+v)];

    // Handle \Top cases
    if (cu[0] == TOP) return cv[0] == TOP ? 0 : 1;
    if (cv[0] == TOP) return -1;

    // Non-\Top case: only the first position where the codes differ matters
    int i = 0;
    while (i < l and cu[i] == cv[i]) i++;
    if (i == l) return 0;
    const int du = (int)(cu[i] >> 1), dv = (int)(cv[i] >> 1);
    const uint32_t bu = cu[i] & 1u, bv = cv[i] & 1u;
    if (du < dv) return bu == 0 ? -1 : 1; // v has epsilon, so if u has 0, -1
    if (du > dv) return bv == 0 ? 1 : -1; // mirror case: u has epsilon, so if v has 0, 1
    return bu < bv ? -1 : 1; // same depth, different bit
}

bool Succinct::eq(int u, int v) const
{
    const uint32_t *cu = &code[(size_t)l*(4+u)];
    const uint32_t *cv = &code[(size_t)l*(4+v)];
    if (cu[0] == TOP or cv[0] == TOP) return cu[0] == cv[0];
    return std::memcmp(cu, cv, sizeof(uint32_t)*l) == 0;
}

void Succinct::stream(std::ostream &out, int target) const
{
    const uint32_t *c = &code[(size_t)l*(4+target)];
    if (c[0] == TOP) {
        out << "\033[1;33mTop\033[m";
    } else {
        out << "[ ";
        int j=0;
        for (int i=0; i<h; i++) {
            if (i>0) out << ",";
            int n=0;
            while (j<l and (int)(c[j]>>1)==i) {
                n++;
                out << (c[j] & 1u);
                j++;
            }
            if (n == 0) out << "ε";
        }
        out << " ]";
    }
}

std::vector<std::byte> Succinct::store(int tgt) const
{
    // per bit: the bit (uint8) and its depth (int); depth -1 at the first bit means Top
    std::vector<std::byte> blob;
    blob.reserve(static_cast<size_t>(l) * (sizeof(uint8_t) + sizeof(int))); // perf: avoid reallocs

    const uint32_t *c = &code[(size_t)l*(4+tgt)];
    for (int i=0; i<l; i++) {
        if (i == 0 and c[0] == TOP) {
            push_le<uint8_t>(blob, 0);
            push_le<int>(blob, -1);
        } else {
            push_le<uint8_t>(blob, (uint8_t)(c[i] & 1u));
            push_le<int>(blob, (int)(c[i] >> 1));
        }
    }

    return blob;
}

void Succinct::load(int tgt, const std::vector<std::byte>& blob)
{
    size_t cur = 0;

    uint32_t *c = &code[(size_t)l*(4+tgt)];
    for (int i = 0; i < l; ++i) {
        const uint32_t b = read_le_advance<uint8_t>(blob, cur) != 0 ? 1u : 0u;
        const int d = read_le_advance<int>(blob, cur);
        c[i] = (i == 0 and d == -1) ? TOP : succ_code(d, b);
    }
}


/**
 * Strahler measures: see gpm.hpp. The successor (prog_tmp) and the comparison follow the
 * STRPM solver by Vanessa Flügel and Guillermo A. Pérez; the truncation (trunc_tmp) is new.
 * A label is stored as one code per position (see gpm.hpp), so that copy and compare are
 * plain loops over integers.
 */

static int
strahler_floor_log2(unsigned long long x)
{
    int y = 0;
    while (x > 1) { x >>= 1; y++; }
    return y;
}

// number of distinct values in the non-decreasing range [first, last)
static inline int
count_distinct_nondecreasing(std::vector<int>::const_iterator first, std::vector<int>::const_iterator last)
{
    if (first == last) return 0;
    int distinct = 1;
    for (auto it = first + 1; it != last; ++it) if (*it != *(it - 1)) ++distinct;
    return distinct;
}

Strahler::Strahler(const Game &game, const int player)
{
    this->player = player;
    this->N = game.nodecount();
    const int max_pr = game.priority(game.nodecount()-1);
    const int h0 = (max_pr/2)+1;
    const int h1 = (max_pr+1)/2;
    const int depth = player == 0 ? h0 : h1;
    this->h = depth + 1; // levels 0..h-2 (as in STRPM)
    this->t = strahler_floor_log2(N);
    this->k = std::min(t + 1, h); // (n,h)-universal (see gpm.hpp); k-1 non-empty strings
                                   // fit in the h-1 levels
    if (depth == 0) { k = 1; t = 0; } // no levels: only the empty label
    this->cap = std::max(0, k - 1 + t);
    this->stride = cap + 1;
    this->END = (uint32_t)h + 1; // above every level (0..h-2)
    this->TOP = 2*END + 1; // above every code of a label
    code.assign((size_t)stride*(N+4), END);
    tb.reserve(cap + 1);
    td.reserve(cap + 1);
    for (int i=-4; i<N; i++) copy_bot(i);
}

void
Strahler::set(const Measures &other)
{
    // (not via operator=: Measures::operator= calls set)
    const Strahler *src = dynamic_cast<const Strahler*>(&other);
    player = src->player;
    k = src->k;
    t = src->t;
    h = src->h;
    cap = src->cap;
    stride = src->stride;
    N = src->N;
    END = src->END;
    TOP = src->TOP;
    code = src->code;
}

void Strahler::copy(int from, int to)
{
    std::memcpy(&code[(size_t)stride*(to+4)], &code[(size_t)stride*(from+4)], sizeof(uint32_t)*stride);
}

void Strahler::copy_from(const Measures& mother, int from, int to)
{
    const Strahler *other = dynamic_cast<const Strahler*>(&mother);
    std::memcpy(&code[(size_t)stride*(to+4)], &other->code[(size_t)stride*(from+4)], sizeof(uint32_t)*stride);
}

void Strahler::copy_bot(int to)
{
    // the leftmost leaf: t+1 zeros at level 0, then a single zero at levels 1..k-2
    uint32_t *c = &code[(size_t)stride*(to+4)];
    for (int i=0; i<cap; i++) c[i] = i <= t ? 0 : (uint32_t)(i - t);
    c[cap] = END;
}

void Strahler::copy_top(int to)
{
    uint32_t *c = &code[(size_t)stride*(to+4)];
    std::fill_n(c, stride, END);
    c[0] = TOP;
}

bool Strahler::is_top(int tgt) const
{
    return code[(size_t)stride*(tgt+4)] == TOP;
}

void Strahler::load_tmp(int slot)
{
    const uint32_t *c = &code[(size_t)stride*(slot+4)];
    tb.resize(cap);
    td.resize(cap);
    int n = 0;
    for (; c[n] != END; n++) { // (the last code of a slot is END)
        const bool one = c[n] > END;
        tb[n] = one ? 1 : 0;
        td[n] = (int)(one ? 2*END - c[n] : c[n]);
    }
    tb.resize(n);
    td.resize(n);
}

void Strahler::store_tmp(int slot)
{
    if (!td.empty() and td[0] == -1) {
        copy_top(slot);
        return;
    }
    uint32_t *c = &code[(size_t)stride*(slot+4)];
    const int n = (int)td.size();
    for (int i=0; i<n; i++) {
        assert(td[i] >= 0 and (uint32_t)td[i] < END);
        c[i] = tb[i] ? 2*END - (uint32_t)td[i] : (uint32_t)td[i];
    }
    std::fill(c + n, c + stride, END);
}

int Strahler::skip_until_next_level(int i)
{
    while (((long)i == (long)td.size() - 1) || (i >= 0 && td[i] == td[i+1])) {
        tb[i] = 0;
        i--;
    }
    return i;
}

void Strahler::prog_tmp(int pindex)
{
    if (td.empty()) { // the only label (no levels)
        tb.assign(1, 0);
        td.assign(1, -1);
        return;
    }
    if (td[0] == -1) return; // already Top

    bool skipLevel = false;
    int i = (int)td.size() - 1;
    // skip bits "below p"
    while (i >= 0 && td[i] > pindex) {
        tb[i] = 0;
        i--;
    }
    // number of non-empty strings (NES) up to here
    int nes = count_distinct_nondecreasing(td.begin(), td.begin() + i + 1);
    if (nes - 1 == k - 1) {
        // no next sibling on this layer
        i = skip_until_next_level(i);
    }

    while (i >= 0) {
        // check if there was a level change
        if ((size_t)i == td.size()-1 || td[i] != td[i+1]) {
            // non-leading bits (NLB): the length minus one leading bit per non-empty string
            int nlb = (i + 1) - nes;
            nes--;
            if (nlb < t) {
                int new_index = std::min(((size_t)(i+1) == td.size()) ? td[i] : td[i+1] - 1, pindex);
                i++;
                // either add one at the end of the existing string (one more NLB), or a new string
                if (i != 0 && ((size_t)i == td.size() || td[i-1] == new_index)) nlb++;
                else nes++;

                int total_nes = std::min((nes + h - 1 - new_index), k - 1);
                if ((long)(total_nes + t) > (long)tb.size()) {
                    size_t newBits = total_nes + t - td.size();
                    tb.insert(tb.end(), newBits, 0);
                    tb[i] = 1;
                    for (size_t j = i; j < td.size(); j++) td[j] = new_index;
                    td.insert(td.end(), newBits, new_index);
                    i += t - nlb + 1;
                } else {
                    tb[i] = 1;
                    td[i] = new_index;
                    int j = 1;
                    while (nlb + j <= t) {
                        tb[i + j] = 0;
                        td[i + j] = new_index;
                        j++;
                    }
                    i += j; // the last changed position
                }
                break;
            }
            // the current string is only a leading bit (all NLB are in the strings before it)
            else if (i == 0 || td[i-1] != td[i]) {
                tb[i] = 0;
                i--;
                continue;
            }
        }
        // here NLB == t
        if (tb[i] == 0) {
            if (i == 0 || td[i - 1] != td[i]) {
                // the 0 is the first bit, or the first bit of its level
                int strings_after_current = count_distinct_nondecreasing(td.begin() + i, td.end());
                if (strings_after_current == (h-1) - td[i]) {
                    // all strings after the current level are non-empty: no sibling on this layer
                    i--;
                    continue;
                } else {
                    // the string is 01^j, so the new string would be empty: one NES less
                    nes--;
                    skipLevel = true;
                }
                break;
            } else {
                break;
            }
        } else {
            tb[i] = 0;
            i--;
        }
    }

    // reached the root: Top
    if (i == -1) {
        if (td[0] == 0) {
            td[0] = -1;
            return;
        } else {
            // there are still empty strings "left" that can be filled instead
            i = 0;
            tb[i] = 1;
            td[i] = std::min(td[i], pindex+1) - 2;
            if ((long)(t + k - 1) > (long)tb.size()) {
                size_t newBits = t + k - 1 - td.size();
                tb.insert(tb.end(), newBits, 0);
                for (size_t j = i; j < td.size(); j++) td[j] = td[i];
                td.insert(td.end(), newBits, td[i]);
            }
            skipLevel = true;
            nes--;
        }
    }

    // we can only fill as many NES as we have height
    int available_nes = h - 1 - (td[(skipLevel ? i : i-1)] + 1);
    int no_of_needed_nes = std::min((k-1) - (nes+1), available_nes);

    if (no_of_needed_nes == 0) {
        // remove everything after the current position
        td.resize(i);
        tb.resize(i);
    } else {
        int set_index = td[(skipLevel ? i : i-1)] + 1;
        // fill up with the next level as long as we still have bits
        while ((long)tb.size() >= (long)(i + no_of_needed_nes)) {
            td[i] = set_index;
            i++;
        }
        // then the rest of the bits one level each
        while ((long)i < (long)tb.size()) {
            set_index++;
            td[i] = set_index;
            i++;
        }
    }
}

void Strahler::trunc_tmp(int pindex)
{
    if (!td.empty() and td[0] == -1) return; // Top
    if (k == 1) {
        // a tree with k = 1 has a single leaf, the bottom label (t zeros at level 0)
        tb.assign(cap, 0);
        td.assign(cap, 0);
        return;
    }
    // keep the strings at levels <= pindex
    int i = (int)td.size() - 1;
    while (i >= 0 && td[i] > pindex) i--;
    tb.resize(i+1);
    td.resize(i+1);
    // then the smallest continuation (0x < epsilon): all remaining non-leading bits as zeros
    // in one string at level pindex+1, then a single zero at each next level, as long as
    // there are non-empty strings and levels left
    const int nes = count_distinct_nondecreasing(td.begin(), td.end());
    const int nlb = (int)td.size() - nes;
    const int more = std::min((k-1) - nes, (h-2) - pindex);
    for (int j = 0; j < more; j++) {
        const int n = j == 0 ? 1 + t - nlb : 1;
        for (int x = 0; x < n; x++) {
            tb.push_back(0);
            td.push_back(pindex + 1 + j);
        }
    }
}

void Strahler::see(int tgt, int priority)
{
    if (is_top(tgt)) return;
    const int pindex = player == 0 ? (h-1)-(priority+1)/2-1 : (h-1)-priority/2-1;
    load_tmp(tgt);
    if ((priority&1) == player) prog_tmp(pindex);
    else trunc_tmp(pindex);
    store_tmp(tgt);
}

void Strahler::inc(int tgt)
{
    see(tgt, player);
}

/**
 * Return -1 if value[u] < value[v], 0 if equal, 1 if value[u] > value[v]
 */
/**
 * At the first position where two labels differ: if one has bit 0 at level d and the other
 * ended its string at level d (so it has a higher level, or nothing, at that position), then
 * the first is smaller (0x < epsilon), and with bit 1 it is larger (epsilon < 1x); with the
 * same level the bits decide. The codes d < (higher levels and END) < 2*END-d give exactly
 * this, so the labels compare as their codes do.
 */
int Strahler::compare(int u, int v) const
{
    const uint32_t *cu = &code[(size_t)stride*(u+4)], *cv = &code[(size_t)stride*(v+4)];
    for (int i=0; i<stride; i++) {
        if (cu[i] != cv[i]) return cu[i] < cv[i] ? -1 : 1;
    }
    return 0;
}

bool Strahler::eq(int u, int v) const
{
    return std::memcmp(&code[(size_t)stride*(u+4)], &code[(size_t)stride*(v+4)], sizeof(uint32_t)*stride) == 0;
}

void Strahler::stream(std::ostream &out, int tgt) const
{
    if (is_top(tgt)) {
        out << "\033[1;33mTop\033[m";
        return;
    }
    const uint32_t *c = &code[(size_t)stride*(tgt+4)];
    out << "[ ";
    int j = 0;
    for (int level=0; level<h-1; level++) {
        if (level > 0) out << ",";
        int n = 0;
        for (; c[j] != END; j++, n++) {
            const bool one = c[j] > END;
            if ((int)(one ? 2*END - c[j] : c[j]) != level) break;
            out << (one ? 1 : 0);
        }
        if (n == 0) out << "ε";
    }
    out << " ]";
}

std::vector<std::byte> Strahler::store(int tgt) const
{
    // the label length (-1 for Top), then per position its bit and level
    std::vector<std::byte> blob;
    const uint32_t *c = &code[(size_t)stride*(tgt+4)];
    if (c[0] == TOP) {
        push_le<int>(blob, -1);
        return blob;
    }
    int n = 0;
    while (c[n] != END) n++;
    push_le<int>(blob, n);
    for (int i=0; i<n; i++) {
        const bool one = c[i] > END;
        push_le<uint8_t>(blob, one ? 1 : 0);
        push_le<int>(blob, (int)(one ? 2*END - c[i] : c[i]));
    }
    return blob;
}

void Strahler::load(int tgt, const std::vector<std::byte>& blob)
{
    size_t cur = 0;
    const int n = read_le_advance<int>(blob, cur);
    if (n == -1) {
        copy_top(tgt);
        return;
    }
    uint32_t *c = &code[(size_t)stride*(tgt+4)];
    for (int i=0; i<n; i++) {
        const uint8_t b = read_le_advance<uint8_t>(blob, cur);
        const int d = read_le_advance<int>(blob, cur);
        c[i] = b ? 2*END - (uint32_t)d : (uint32_t)d;
    }
    std::fill(c + n, c + stride, END);
}


GPMSolver::GPMSolver(Oink& oink, Game& game, MeasureKind kind) : Solver(oink, game), measure_kind(kind)
{
}

GPMSolver::~GPMSolver()
{
}

bool
GPMSolver::lift(Measures &pm, const int player, int v) {
    bool updated = false;
    if (owner(v) == player) {
        // Take the highest successor
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            int to = *curedge;
            if (!G[to]) continue;
            pm.copy(to, -1);
            pm.see(-1, priority(v));
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "Lift option (takes highest) " << label_vertex(v) << " to " << label_vertex(to) << " ";
                pm.stream(logger, -1);
                logger << std::endl;
            }
#endif
            if ((!updated || pm.compare(-1, -4) > 0)) {
                pm.copy(-1, -4);
                updated = true;
            }
        }
    } else {
        // Take the lowest successor
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            int to = *curedge;
            if (!G[to]) continue;
            pm.copy(to, -1);
            pm.see(-1, priority(v));
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "Lift option (takes lowest) " << label_vertex(v) << " to " << label_vertex(to) << " ";
                pm.stream(logger, -1);
                logger << std::endl;
            }
#endif
            if ((!updated || pm.compare(-1, -4) <= 0)) {
                pm.copy(-1, -4);
                updated = true;
            }
        }
    }
    return updated;
}

bool
GPMSolver::update(Measures &pm, const int player) {
    auto &Q = player == 0 ? Q0 : Q1;
    //auto &R = player == 0 ? R0 : R1;
    auto &T = player == 0 ? T0 : T1;
    //auto &B = player == 0 ? B0 : B1;
    auto &L = player == 0 ? L0 : L1;

    QQ.reset();
    L.reset();
    for (auto v = Q.find_first(); v != bitset::npos; v = Q.find_next(v)) {
        lift_attempts++;
        if (lift(pm, player, v)) {
            auto res = pm.compare(-4, v);
            if (res < 0) {
#ifndef NDEBUG
                if (trace) {
                    logger << "\033[1m" << player << "-lifting\033[m \033[36m" << label_vertex(v) << "\033[m from \033[1m";
                    pm.stream(logger, v);
                    logger << " to ";
                    pm.stream(logger, -4);
                    logger << std::endl;
                }
#endif
                //LOGIC_ERROR;
            }
            if (res > 0) {
                lifts++;
#ifndef NDEBUG
                if (trace) {
                    logger << "\033[1m" << player << "-lifting\033[m \033[36m" << label_vertex(v) << "\033[m from \033[1m";
                    pm.stream(logger, v);
                    logger << "\033[m to \033[1m";
                    pm.stream(logger, -4);
                    logger << "\033[m" << std::endl;
                }
#endif
                L[v] = 1;
                // actually lift
                pm.copy(-4, v);
                // update B and R
                if (pm.is_top(-4)) {
                    T[v] = true;
                    //R[v] = false;
                }
                //B[v] = false;
                // update QQ
                for (auto curedge = ins(v); *curedge != -1; curedge++) {
                    int from = *curedge;
                    if (!T[from]) QQ[from] = true;
                }
            }
        }
    }

    Q.swap(QQ);
    return Q.any();
}


/**
 * Attract as player <pl> via <v> to <Z>, vertices in <R>. Only escape to <R>.
 * Add attracted vertices to <Z> and to queue <Q>.
 */
void
GPMSolver::attractVertices(int pl, int v, bitset &R, bitset &Z)
{
    // attract vertices with an edge to <v>
    for (auto curedge = ins(v); *curedge != -1; curedge++) {
        int from = *curedge;
        if (!Z[from] && R[from]) {
            if (owner(from) != pl) {
                // check if opponent can escape
                bool escapes = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (R[to] && !Z[to]) {
                        escapes = true;
                        break;
                    }
                }
                if (escapes) continue;
            }
            // attract
            Z[from] = true;
            Queue.push(from);
        }
    }
}


void
GPMSolver::shortcuts(const int player, Measures &pm0, Measures &pm1)
{
    auto& Q = player == 0 ? Q0 : Q1;
    auto& Qo = player == 0 ? Q1 : Q0;
    auto& T = player == 0 ? T0 : T1;
    auto& To = player == 0 ? T1 : T0;
    auto& pm = player == 0 ? pm0 : pm1;
    auto& pm_o = player == 0 ? pm1 : pm0;
    auto& L = player == 0 ? L0 : L1;

/*
    // Attract towards Bottom!
    QQ = B;
    for (auto v = QQ.find_first(); v != bitset::npos; v = QQ.find_next(v)) {
        logger << "bottom: " << label_vertex(v) << std::endl;
        Queue.push(v);
    }
    while (Queue.nonempty()) {
        const int v = Queue.pop();
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            int from = *curedge;
            if (QQ[from] || !R[from]) continue; // don't attract
            if (owner(from) == player) {
                // see we can escape (not have to play to bottom)
                bool escapes = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (!R[to] || QQ[to]) continue; // not a potential escape
                    // it is only a valid escape if the measure obtained by going from
                    // <from> to <to> actually is >= the current measure...
                    pm.copy(to, -1);
                    pm.see(-1, priority(from));
                    if (pm.compare(-1, from) >= 0) {
                        escapes = true;
                        break;
                    }
                }
                if (escapes) continue;
            }
            // attract
            QQ[from] = true;
            Queue.push(from);
        }
    }
    QQ ^= G;
    QQ -= T;
    for (auto v = QQ.find_first(); v != bitset::npos; v = QQ.find_next(v)) {
        logger << "player " << player << " should win " << label_vertex(v) << std::endl;
    }
*/




    // logger<<"running shortcuts for player " << player << std::endl; 

    // now find vertices that are "stable"
    // any vertex in L is lifted previously, so it's increasing for L
    // also add all in T
    L |= T;
    for (auto v = L.find_first(); v != bitset::npos; v = L.find_next(v)) {
        Queue.push(v);
    }
    while (Queue.nonempty()) {
        const int v = Queue.pop();
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            int from = *curedge;
            if (L[from] || !G[from]) continue; // don't attract
            if (owner(from) != player) {
                // see if opponent can escape
                bool escapes = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (!G[to] || L[to]) continue; // not a potential escape
                    // it is only a valid escape if the measure obtained by going from
                    // <from> to <to> actually is == the current measure...
                    pm.copy(to, -1);
                    pm.see(-1, priority(from));
                    if (pm.compare(-1, from) <= 0) {
                        escapes = true;
                        break;
                    }
                }
                if (escapes) continue;
            }
            // attract
            L[from] = true;
            Queue.push(from);
        }

    }
    L ^= G;
    for (auto v = L.find_first(); v != bitset::npos; v = L.find_next(v)) {
        // logger << "player " << (1-player) << " should win " << label_vertex(v) << std::endl;
        pm_o.copy_top(v);
        To[v] = true;
        //Ro[v] = false;
        Qo[v] = false;
        //Bo[v] = false;
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            int from = *curedge;
            if (G[from] && !To[from]) Qo[from] = true;
        }
    }
}


void
GPMSolver::solve(Measures &pm, const int player)
{
    // extract strategy
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (!pm.is_top(v)) {
            // won by ~X
            str[v] = -1;
            if (owner(v) == player) {
                // find lowest prog...
                for (auto curedge = outs(v); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (G[to] && !pm.is_top(to)) {
                        if (str[v] == -1) {
                            pm.copy(to, -1);
                            str[v] = to;
                        } else if (pm.compare(to, -1) < 0) {
                            pm.copy(to, -1);
                            str[v] = to;
                        }
                    }
                }

                if (str[v] == -1) LOGIC_ERROR;
            }
        }
    }
    // mark as solved
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (!pm.is_top(v)) {
#ifndef NDEBUG
            if(trace) {
                logger << "player " << player << " wins " << label_vertex(v) << " with strategy " << label_vertex(str[v]) << " and ";
                pm.stream(logger, v);
                logger << std::endl;
            }
#endif
            Solver::solve(v, player, str[v]);
            G[v] = false;
        }
    }
    Solver::flush();
}

bool
GPMSolver::parseOptions(std::string& opts)
{
    if (opts.empty()) return true; // keep the defaults (or the measure kind of gpm-small etc.)
    std::string remaining = opts;

    do {
        size_t comma_index = remaining.find(',');
        std::string opt = remaining.substr(0, comma_index);
        if (opt == "noshortcuts" || opt == "no-shortcuts") {
            use_shortcuts = false;
        } else {
            std::optional<MeasureKind> parsed = parse_measure_kind(opt);
            if (parsed.has_value()) {
                measure_kind = *parsed;
            } else {
                logger << "unrecognised option for GPM solver, possible options: ";
                stream_measure_kinds(logger);
                logger << ", no-shortcuts (separated by commas and no spaces)" << std::endl;
                return false;
            }
        }
        if (comma_index == std::string::npos) {
            remaining = remaining.substr(remaining.size());
        } else {
            remaining = remaining.substr(comma_index + 1);
        }
    } while (!remaining.empty());

    if (trace) {
        logger << "selected ";
        stream_measure_kind(logger, measure_kind);
        logger << " for GPM" << std::endl;
        if (!use_shortcuts) {
            logger << "disabled shortcuts for GPM" << std::endl;
        }
    }

    return true;
}

void
GPMSolver::run()
{
    str = new int[nodecount()];
    G = disabled;
    G.flip();

    //Q = disabled;
    // B0 = B1 = G;
    // R0 = R1 = G;

    Q0 = Q1 = G;
    QQ.resize(nodecount());
    Queue.resize(nodecount());
    T0.resize(nodecount());
    T1.resize(nodecount());
    L0.resize(nodecount());
    L1.resize(nodecount());

    Measures *pm0 = new_measure(measure_kind, game, 0);
    Measures *pm1 = new_measure(measure_kind, game, 1);

    bool c0 = true;
    bool c1 = true;

    long iterations = 0;

    while (true) {
        iterations++;
        if (c0) {
            c0 = update(*pm0, 0);
            if (use_shortcuts) {
                if (!c0 || (iterations%10) == 0) shortcuts(0, *pm0, *pm1);
            }
        }

        if (c1) {
            c1 = update(*pm1, 1);
            if (use_shortcuts) {
                if (!c1 || (iterations%10) == 0) shortcuts(1, *pm0, *pm1);
            }
        }

        if (!c0 && !c1) break;
    }

    solve(*pm1, 0);
    solve(*pm0, 1);

    logger << "solved with " << lifts << " lifts, " << lift_attempts << " lift attempts, " << iterations << " iterations." << std::endl;

    delete[] str;
    delete pm0;
    delete pm1;
}

std::optional<MeasureKind> parse_measure_kind(std::string& opts) {
    if (opts.empty()) {
        return MeasureKind::Ordered;
    } else if (opts == "opm" || opts == "qpt" || opts == "ordered" || opts == "OPM" || opts == "QPT" || opts == "ORDERED") {
        return { MeasureKind::Ordered };
    } else if (opts == "spm" || opts == "small" || opts == "SPM" || opts == "SMALL") {
        return { MeasureKind::Small };
    } else if (opts == "sspm"  || opts == "succinct" || opts == "SSPM" || opts == "SUCCINCT") {
        return { MeasureKind::Succinct };
    } else if (opts == "strpm" || opts == "strahler" || opts == "STRPM" || opts == "STRAHLER") {
        return { MeasureKind::Strahler };
    }

    return std::nullopt;
}

Measures *new_measure(const MeasureKind kind, const Game &game, const int player) {
    switch (kind) {
        case MeasureKind::Small:
            return new Small(game, player);
        case MeasureKind::Ordered:
            return new Ordered(game, player);
        case MeasureKind::Succinct:
            return new Succinct(game, player);
        case MeasureKind::Strahler:
            return new Strahler(game, player);
        default:
            THROW_ERROR("Measure kind exists but is not implemented");
    }
}

void stream_measure_kind(std::ostream& logger, MeasureKind kind) {
    switch (kind) {
        case MeasureKind::Small:
            logger << "Small Progress Measures (SPM)";
            break;
        case MeasureKind::Ordered:
            logger << "Ordered Progress Measures (OPM)";
            break;
        case MeasureKind::Succinct:
            logger << "Succinct Progress Measures (SSPM)";
            break;
        case MeasureKind::Strahler:
            logger << "Strahler Progress Measures (STRPM)";
            break;
        default:
            THROW_ERROR("Measure kind exists but is not implemented");
            break;
    }
}

void stream_measure_kinds(std::ostream& logger) {
    logger << "small, ordered, succinct, strahler";
}

}
