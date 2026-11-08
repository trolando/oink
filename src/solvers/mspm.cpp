/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
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
#include <cmath>
#include <string>

#include "mspm.hpp"

namespace pg {

MSPMSolver::MSPMSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

MSPMSolver::~MSPMSolver()
{
}

/**
 * Whether measure <a> is less than <b> up to and including priority <d>.
 */
bool
MSPMSolver::pm_less(const int *a, const int *b, int d) const
{
    if (b[P] == -1) return a[P] != -1;
    if (a[P] == -1) return false;
    const int start = ((k-1)&1) == P ? k-1 : k-2; // the highest position of parity <P>
    for (int i=start; i>=d; i-=2) {
        if (a[i] != b[i]) return a[i] < b[i];
    }
    return false;
}

/**
 * The least measure for a vertex with priority <d> and a successor with measure <src>:
 * at least <src> up to priority <d>, and higher than <src> if <d> has parity <P>.
 */
void
MSPMSolver::Prog(int *dst, const int *src, int d) const
{
    if (src[P] == -1) {
        dst[P] = -1;
        return;
    }
    int i = P;
    for (; i<d; i+=2) dst[i] = 0;
    int carry = (i == d) ? 1 : 0; // only increase for a priority of parity <P>
    for (; i<k; i+=2) {
        const int v = src[i] + carry;
        if (v > counts[i]) {
            dst[i] = 0;
            carry = 1;
        } else {
            dst[i] = v;
            carry = 0;
        }
    }
    if (carry) dst[P] = -1;
}

/**
 * Lift <v>: the highest Prog over its successors if <v> is of <P>, otherwise the lowest.
 */
bool
MSPMSolver::lift(int v)
{
    int *pm = &pms[(size_t)k*v];
    if (pm[P] == -1) return false;
    lift_attempt++;

    const int d = priority(v);
    const bool take_max = owner(v) == P;
    bool have = false;
    for (auto curedge = outs(v); *curedge != -1; curedge++) {
        const int w = *curedge;
        if (depth[w] == -2) continue; // disabled
        Prog(tmp.data(), &pms[(size_t)k*w], d);
        if (!have or (take_max ? pm_less(best.data(), tmp.data(), d) : pm_less(tmp.data(), best.data(), d))) {
            for (int i=P; i<k; i+=2) best[i] = tmp[i];
            have = true;
        }
    }
    if (!have or !pm_less(pm, best.data(), d)) return false;

    for (int i=P; i<k; i+=2) pm[i] = best[i];
    lift_count++;
    return true;
}

/**
 * Lift the vertices of context <D> (in <W>) until no vertex can be lifted, or until a
 * vertex reaches Top; returns that vertex, or -1.
 */
int
MSPMSolver::lift_context(int D, const std::vector<int> &W)
{
    for (const int w : W) {
        if (depth[w] == D and !is_top(w) and !dirty[w]) {
            dirty[w] = true;
            Q.push(w);
        }
    }
    int found = -1;
    while (Q.nonempty()) {
        const int w = Q.pop();
        dirty[w] = false;
        if (found != -1 or depth[w] != D) continue; // (after Top: just empty the queue)
        if (!lift(w)) continue;
        if (is_top(w)) {
            found = w;
            continue;
        }
        for (auto curedge = ins(w); *curedge != -1; curedge++) {
            const int u = *curedge;
            if (depth[u] == D and !dirty[u] and !is_top(u)) {
                dirty[u] = true;
                Q.push(u);
            }
        }
    }
    return found;
}

/**
 * The attractor for <pl> to <targets> within context <D> (the vertices v with depth[v] == D),
 * only attracting vertices of priority at most <max_prio> (if >= 0) and not in <exclude>.
 * The opponent's vertices are attracted if all their successors in the context are. Marks
 * the attractor in <in> and adds it to <out> (targets included); with <strategy>, sets the
 * strategy of the attracted vertices of <pl>.
 */
void
MSPMSolver::attract(int pl, int D, const std::vector<int> &targets, int max_prio, bitset &in, const bitset *exclude, std::vector<int> &out, bool strategy)
{
    for (const int t : targets) {
        if (in[t]) continue;
        in[t] = true;
        out.push_back(t);
        Q.push(t);
    }
    std::vector<int> touched;
    while (Q.nonempty()) {
        const int n = Q.pop();
        for (auto curedge = ins(n); *curedge != -1; curedge++) {
            const int u = *curedge;
            if (depth[u] != D or in[u]) continue;
            if (exclude != nullptr and (*exclude)[u]) continue;
            if (max_prio >= 0 and priority(u) > max_prio) continue;
            if (owner(u) == pl) {
                if (strategy) str[u] = n;
            } else {
                if (cnt[u] < 0) {
                    int c = 0;
                    for (auto e = outs(u); *e != -1; e++) if (depth[*e] == D) c++;
                    cnt[u] = c;
                    touched.push_back(u);
                }
                if (--cnt[u] > 0) continue;
                if (strategy) str[u] = -1;
            }
            in[u] = true;
            out.push_back(u);
            Q.push(u);
        }
    }
    for (const int u : touched) cnt[u] = -1;
}

/**
 * Options: "even" or "odd" to lift the measure of that player (default: of the player with
 * the smaller measure domain).
 */
bool
MSPMSolver::parseOptions(std::string& opts)
{
    if (opts.empty()) return true;
    if (opts == "even") {
        force_player = 0;
    } else if (opts == "odd") {
        force_player = 1;
    } else {
        logger << "unrecognised option for mspm: " << opts << " (options: even, odd)" << std::endl;
        return false;
    }
    return true;
}

void
MSPMSolver::run()
{
    const int n = nodecount();
    k = priority(n-1)+1;
    if (k < 2) k = 2; // positions 0 and 1 hold the Top marks of the two measures

    pms.assign((size_t)k*n, 0);
    counts.assign(k, 0);
    tmp.assign(k, 0);
    best.assign(k, 0);
    str.assign(n, -1);
    depth.assign(n, 0);
    cnt.assign(n, -1);
    Q.resize(n);
    dirty.resize(n);
    inA.resize(n);
    inB.resize(n);

    std::vector<int> all;
    for (int v=0; v<n; v++) {
        if (disabled[v]) {
            depth[v] = -2;
        } else {
            counts[priority(v)]++;
            all.push_back(v);
        }
    }

    // lift the measure with the smaller domain (as in the paper, either one gives both
    // winning regions and strategies)
    double size[2] = {0, 0};
    for (int i=0; i<k; i++) size[i&1] += std::log((double)counts[i]+1);
    P = force_player >= 0 ? force_player : (size[0] < size[1] ? 0 : 1);
    if (trace) logger << "lifting the measure of player " << P << std::endl;

    std::vector<Frame> stack;
    stack.push_back(Frame{0, std::move(all), {}, {}, false});

    while (!stack.empty()) {
        if (stack.back().child) {
            // the context REM is done: remove the attractor of RES and the Top vertices of REM
            Frame &F = stack.back();
            F.child = false;
            std::vector<int> dom = std::move(F.res);
            F.res.clear();
            for (const int u : F.rem) {
                depth[u] = F.D; // back in this context (also if it was removed in REM)
                if (is_top(u)) dom.push_back(u);
            }
            F.rem.clear();
            std::vector<int> A;
            attract(P, F.D, dom, -1, inA, nullptr, A, true);
            for (const int u : A) {
                pms[(size_t)k*u+P] = -1;
                depth[u] = -1;
                inA[u] = false;
            }
            continue;
        }

        Frame &F = stack.back();
        F.W.erase(std::remove_if(F.W.begin(), F.W.end(), [&](int u) { return depth[u] != F.D; }), F.W.end());
        if (F.W.empty()) {
            stack.pop_back();
            continue;
        }

        const int v = lift_context(F.D, F.W);
        if (v == -1) {
            stack.pop_back(); // no vertex reaches Top: the context is won by the opponent
            continue;
        }

        const int pr = priority(v);
        if (trace) logger << "Top in context " << F.D << ": " << label_vertex(v) << "/" << pr << std::endl;

        if (owner(v) == P) {
            // play to the successor in the context with the highest measure up to <pr>
            int u_best = -1;
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                const int u = *curedge;
                if (depth[u] != F.D) continue;
                if (u_best == -1 or pm_less(&pms[(size_t)k*u_best], &pms[(size_t)k*u], pr)) u_best = u;
            }
            if (u_best == -1) LOGIC_ERROR; // (it reached Top from its successors in the context)
            str[v] = u_best;
        }

        // RES: the attractor of <P> to <v> of the vertices of priority at most <pr>
        std::vector<int> res;
        attract(P, F.D, {v}, pr, inA, nullptr, res, true);
        for (const int u : res) pms[(size_t)k*u+P] = -1;

        // IRR: the attractor of the opponent to the higher priorities
        std::vector<int> high, irr;
        for (const int u : F.W) {
            if (depth[u] == F.D and !inA[u] and priority(u) > pr) high.push_back(u);
        }
        attract(1-P, F.D, high, -1, inB, &inA, irr, false);

        // REM: the rest is the context of the child
        Frame child{F.D+1, {}, {}, {}, false};
        for (const int u : F.W) {
            if (depth[u] != F.D or inA[u] or inB[u]) continue;
            depth[u] = child.D;
            child.W.push_back(u);
        }
        for (const int u : res) inA[u] = false;
        for (const int u : irr) inB[u] = false;
        F.res = std::move(res);
        F.rem = child.W;
        F.child = true;
        stack.push_back(std::move(child)); // (invalidates <F>)
    }

    // the vertices with Top are won by <P> (with the strategy found), the others by the
    // opponent, who plays to a successor with the lowest measure. (Not via <disabled>:
    // Solver::solve adds to it.)
    for (int v=0; v<n; v++) {
        if (depth[v] == -2) continue;
        if (is_top(v)) {
            Solver::solve(v, P, owner(v) == P ? str[v] : -1);
            continue;
        }
        int s = -1;
        if (owner(v) != P) {
            const int d = priority(v);
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                const int w = *curedge;
                if (depth[w] == -2) continue;
                Prog(tmp.data(), &pms[(size_t)k*w], d);
                if (s == -1 or pm_less(tmp.data(), best.data(), d)) {
                    for (int i=P; i<k; i+=2) best[i] = tmp[i];
                    s = w;
                }
            }
        }
        Solver::solve(v, 1-P, s);
    }

    logger << "solved with " << lift_count << " lifts, " << lift_attempt << " lift attempts." << std::endl;
}

}
