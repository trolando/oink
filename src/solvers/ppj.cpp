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
#include <algorithm>
#include <climits>
#include "ppj.hpp"

namespace pg {

// Membership snapshots are contiguous. Invalidated and moved members are
// filtered on the next visit, so attraction writes only the region assignment.
void
PPJSolver::move(int v, int p)
{
    const int old = region[v];
    if (old >= 0 && old != SOLVED && --count[old] == 0) nonempty[old&1][old] = false;
    region[v] = p;
    if (p >= 0 && p != SOLVED && count[p]++ == 0) nonempty[p&1][p] = true;
}

void
PPJSolver::storeMembers(int p)
{
    if (members.size()+work.size() > members.capacity()) {
        std::vector<int> fresh;
        fresh.reserve(std::max(members.capacity(), size_t(nodecount())*2+work.size()));
        for (unsigned q = 0; q < length.size(); ++q) {
            const unsigned b = fresh.size();
            for (unsigned i = begin[q], end = i+length[q]; i < end; ++i) {
                const int v = members[i];
                if (region[v] == int(q)) fresh.push_back(v);
            }
            begin[q] = b;
            length[q] = fresh.size()-b;
        }
        members.swap(fresh);
    }
    begin[p] = members.size();
    for (unsigned i = 0; i < work.size(); ++i) {
        const int v = work[i];
        // Native heads are seeded directly from the priority-sorted game.
        if (native[v] != p) members.push_back(v);
    }
    length[p] = members.size()-begin[p];
}

void
PPJSolver::invalidate(int v)
{
    const int p = region[v];
    dirty[native[v]] = true;
    str[v] = -1;
    if (p == native[v]) return; // the head becomes a parameter of the same winner
    move(v, -1);
    invalid.push(v);
    resets++;
}

void
PPJSolver::invalidateClosure(int p, int pl)
{
    while (invalid.nonempty()) {
        const int v = invalid.pop();
        for (const int *e = ins(v); *e != -1; ++e) {
            const int u = *e, r = region[u];
            if (r < 0 || r >= p || (pl >= 0 && (r&1) == pl)) continue;
            if (str[u] == -1 || str[u] == v) invalidate(u);
        }
    }
}

int
PPJSolver::findEscape(int v, int i, int p, int pl)
{
    const int *edges = outs(v);
    for (int to; (to = edges[i]) != -1; ++i) {
        const int r = region[to];
        if (r < p || (r != SOLVED && (r&1) != pl)) return i;
    }
    return -1;
}

void
PPJSolver::grow(int p)
{
    const int pl = p&1;
    work.clear();
    flipped.clear();
    for (int v = headfirst[p]; v < headend[p]; ++v) {
        if (region[v] == -1) move(v, p);
        if (region[v] == p) work.push(v);
    }
    for (unsigned i = begin[p], end = i+length[p]; i < end; ++i) {
        const int v = members[i];
        if (region[v] == p) work.push(v);
    }
    if (work.empty()) return;

    const bool higher = nonempty[pl].find_next(p) != bitset::npos;
    int changed_lo = p;
    for (unsigned i = 0; i < work.size(); ++i) {
        const int v = work[i];
        for (const int *e = ins(v); *e != -1; ++e) {
            const int u = *e, old = region[u];
            if (old >= p) continue;
            int s = -1;
            if (owner(u) == pl) {
                s = v;
            } else {
                if (cnt[u] == -1) {
                    touched.push(u);
                    if (higher) {
                        witness[u] = findEscape(u, 0, p, pl);
                        cnt[u] = witness[u] == -1 ? 0 : 1;
                    } else {
                        cnt[u] = degree[u]-1;
                    }
                } else if (!higher) {
                    --cnt[u];
                } else if (region[outs(u)[witness[u]]] == p) {
                    witness[u] = findEscape(u, witness[u]+1, p, pl);
                    if (witness[u] == -1) cnt[u] = 0;
                }
                if (cnt[u] != 0) continue;
            }
            if (old >= 0 && (old&1) != pl) flipped.push(u);
            const int lo = old >= 0 && (old&1) == pl ? old+1 : 0;
            changed_lo = std::min(changed_lo, lo);
            move(u, p);
            str[u] = s;
            work.push(u);
        }
    }
    // Only changes of winner break a lower justification. Delay these resets
    // until the attractor is complete, so already captured vertices survive.
    for (unsigned i = 0; i < flipped.size(); ++i) invalid.push(flipped[i]);
    invalidateClosure(p, pl);

    // Changes of level make new targets available in intervening regions.
    for (int side = 0; side < 2; ++side) {
        auto q = changed_lo == 0 ? nonempty[side].find_first() : nonempty[side].find_next(changed_lo-1);
        for (; q != bitset::npos && int(q) < p; q = nonempty[side].find_next(q)) dirty[q] = true;
    }

    // Only native heads can be open. Every other member retains a direct
    // justification, so an open region needs no scan of its interior.
    bool closed = true;
    int target = SOLVED;
    for (int v = headfirst[p]; v < headend[p]; ++v) {
        if (region[v] != p) continue;
        if (owner(v) == pl) {
            if (str[v] == -1) {
                for (const int *e = outs(v); *e != -1; ++e) {
                    const int r = region[*e];
                    if (r != SOLVED && r >= p && (r&1) == pl) {
                        str[v] = *e;
                        break;
                    }
                }
            }
            if (str[v] == -1) {
                closed = false;
                break;
            }
        } else {
            for (const int *e = outs(v); *e != -1; ++e) {
                const int r = region[*e];
                if (r != SOLVED && (r < p || (r&1) != pl)) {
                    closed = false;
                    break;
                }
            }
            if (!closed) break;
        }
    }
    if (closed && higher) {
        // External justification edges determine the first level at which
        // the entire closed region can be reused.
        for (unsigned i = 0; i < work.size(); ++i) {
            const int v = work[i];
            if (owner(v) == pl) {
                const int r = region[str[v]];
                if (r > p) target = std::min(target, r);
            } else {
                for (const int *e = outs(v); *e != -1; ++e) {
                    const int r = region[*e];
                    if (r > p) target = std::min(target, r);
                }
            }
        }
    }
    const bool attract = higher || nonempty[1-pl].find_next(p) != bitset::npos;
    if (!closed || target != SOLVED || attract) {
        for (unsigned i = 0; i < touched.size(); ++i) cnt[touched[i]] = -1;
        touched.clear();
    }
    if (!closed) {
        storeMembers(p);
        return;
    }

    if (target == SOLVED) {
        solveDominion(p, attract);
    } else {
        if (trace) logger << "promoting " << p << " to " << target << " (" << work.size() << " vertices)" << std::endl;
        const unsigned promoted = work.size();
        // Filter the target snapshot before moving: a stale entry can refer
        // to a member currently being promoted, and must not be added twice.
        for (unsigned i = begin[target], end = i+length[target]; i < end; ++i) {
            const int v = members[i];
            if (region[v] == target) work.push(v);
        }
        for (unsigned i = 0; i < promoted; ++i) move(work[i], target);
        length[p] = 0;
        storeMembers(target);
        dirty[target] = true;
        for (int side = 0; side < 2; ++side) {
            for (auto q = nonempty[side].find_first(); q != bitset::npos && int(q) < target; q = nonempty[side].find_next(q)) {
                dirty[q] = true;
            }
        }
        promotions++;
    }
}

void
PPJSolver::solveDominion(int p, bool attract)
{
    const int pl = p&1;
    dominions++;
    length[p] = 0;
    if (trace) logger << "dominion " << p << " (" << work.size() << " vertices before attraction)" << std::endl;
    // Solve each vertex when it is attracted, and attract in stack order, as Oink::flush
    // does: the vertex attracted last is expanded next, while its data is still in the cache.
    for (unsigned i = 0; i < work.size(); ++i) {
        const int v = work[i];
        G[v] = false;
        move(v, SOLVED);
        solve(v, pl, str[v]);
    }
    while (attract && work.nonempty()) {
        const int v = work.pop();
        for (const int *e = ins(v); *e != -1; ++e) {
            const int u = *e;
            if (!G[u]) continue;
            int s = v;
            if (owner(u) != pl) {
                if (cnt[u] < 0) {
                    cnt[u] = degree[u];
                    touched.push(u);
                }
                if (--cnt[u] != 0) continue;
                s = -1;
            }
            G[u] = false;
            move(u, SOLVED);
            solve(u, pl, s);
            work.push(u);
        }
    }
    work.clear();
    for (unsigned i = 0; i < touched.size(); ++i) {
        const int u = touched[i];
        degree[u] = cnt[u];
        cnt[u] = -1;
        if (!G[u] || region[u] == -1) continue;
        // (justifications never point to vertices solved in an earlier dominion)
        if (str[u] != -1 && !G[str[u]]) invalidate(u);
    }
    touched.clear();
    invalidateClosure(SOLVED, -1);
    dirty |= nonempty[0];
    dirty |= nonempty[1];
}

void
PPJSolver::run()
{
    const int n = nodecount();
    G = disabled;
    G.flip();
    if (G.none()) return;
    const int maxp = priority(G.find_last());
    dirty.resize(maxp+1);
    nonempty[0].resize(maxp+1);
    nonempty[1].resize(maxp+1);
    count.assign(maxp+1, 0);
    begin.assign(maxp+1, 0);
    length.assign(maxp+1, 0);
    members.reserve(size_t(n)*2);
    headfirst.assign(maxp+1, 0);
    headend.assign(maxp+1, 0);
    // These arrays are initialized on first use. A shared allocation avoids
    // clearing vertex-sized buffers whose entries are immediately overwritten.
    workspace.reset(new int[size_t(n)*6]);
    region = workspace.get();
    str = region+n;
    native = str+n;
    degree = native+n;
    cnt = degree+n;
    witness = cnt+n;
    work.resize(n);
    invalid.resize(n);
    flipped.resize(n);
    touched.resize(n);
    const bool all = disabled.none();
    int top = -1;
    for (int v = n-1; v >= 0; --v) {
        region[v] = SOLVED;
        if (!G[v]) continue;
        str[v] = -1;
        cnt[v] = -1;
        degree[v] = 0;
        const int pr = priority(v);
        if (top == -1 || ((top^pr)&1)) {
            top = pr;
            headend[top] = v+1;
            dirty[top] = true;
        }
        native[v] = top;
        headfirst[top] = v;
        region[v] = -1;
        if (all) {
            degree[v] = game.outcount(v);
        } else {
            for (const int *e = outs(v); *e != -1; ++e) {
                if (G[*e]) degree[v]++;
            }
        }
    }
    while (dirty.any()) {
        const int p = dirty.find_last();
        dirty[p] = false;
        ++iterations;
        grow(p);
    }
    if (G.any()) THROW_ERROR("PPJ search was incomplete");
    logger << "found " << dominions << " dominions, " << promotions << " promotions and " << resets << " resets." << std::endl;
    logger << "solved in " << iterations << " iterations." << std::endl;
}
}
