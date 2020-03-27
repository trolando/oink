/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
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
#include <cassert>
#include <cstring> // for memset

#include "rtl.hpp"

namespace pg {

// TODO
// store all edges in a vector
// in-edges consecutive
// store per vertex: i_first (i_last from next)
// for out: extra indirection?

RTLSolver::RTLSolver(Oink *oink, Game *game) : Solver(oink, game)
{
}

RTLSolver::~RTLSolver()
{
}

/**
 * Attract to <v> in the subgame <R> as player <pl>.
 * Update <Z> and <str>.
 * Adds each attracted vertex to the Todo queue.
 */
void
RTLSolver::attractVertices(const int pr, const int pl, int v, bitset &R, int *str, bitset &Z)
{
    // attract vertices
    for (auto curedge = ins(v); *curedge != -1; curedge++) {
        int from = *curedge;
        if (Z[from]) {
            // already in Z, maybe set strategy
            if (owner(from) == pl and str[from] == -1) str[from] = v;
        } else if (R[from]) {
            // in subgame but not (yet) in Z
            if (owner(from) != pl) {
                // check each exit
                bool escapes = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (R[to] and !Z[to]) {
                        escapes = true;
                        break;
                    }
                }
                if (escapes) continue;
            }
            // attract
            Z[from] = true;
            str[from] = owner(from) == pl ? v : -1;
            Q.push(from);
#ifndef NDEBUG
            // maybe report event
            if (trace >= 3) {
                logger << "\033[1;37mattracted \033[36m" << label_vertex(from) << "\033[m to \033[1;36m" << pr << "\033[m";
                if (owner(from) == pl) logger << " (via " << label_vertex(v) << ")" << std::endl;
                else logger << " (forced)" << std::endl;
            }
#endif
        }
    }
    (void)pr; // suppress warning if NDEBUG is defined
}

/**
 * Attract vertices that are in tangles.
 * Current attracting set is <Z>.
 * Current subgame is <Z> + <R>.
 * Write strategy to <str>.
 * Current attracting vertex is <v>.
 * Attracting for player <pl>.
 * (for trace) Attracting to region with priority <pr>.
 */
bool
RTLSolver::attractTangle(const int t, const int pl, bitset &R, int *str, bitset &Z)
{
    /**
     * Check if tangle in T_alpha, i.e., our parity, not deleted.
     */

    const int tangle_pr = tpr[t];
    if (tangle_pr == -1) return false; // deleted tangle
    if ((tangle_pr&1) != pl) return false; // not of our parity

    /**
     * Check if tangle is contained in Z+R.
     * We require at least one vertex in Z\R.
     */

    {
        bool can_attract = false;
        int *ptr, x;
        ptr = tv[t];
        while ((x=*ptr++) != -1) {
            ptr++;
            if (disabled[x]) {
                // on-the-fly detect disabled tangles
                // now mark permanently as disabled and break
                tpr[t] = -1;
                return false;
            } else if (Z[x]) {
                continue; // already attracted
            } else if (R[x]) {
                can_attract = true; // has vertices not yet attracted
            } else {
                return false; // not contained in Z+R
            }
        }
        if (!can_attract) return false; // either no vertices in R\Z or strategy leaves R+Z
    }

    /**
     * Check if the tangle can escape to R\Z.
     * Future TODO: actually use the edges instead... (they can be disabled)
     */

    {
        int *ptr, x;
        ptr = tout[t];
        while ((x=*ptr++) != -1) {
            if (R[x] and !Z[x]) return false;
        }
    }

    /**
     * Attract!
     * Future TODO: only add vertices to Q that have ingoing edges (precompute)
     */

    {
        int *ptr, x;
        ptr = tv[t];
        while ((x=*ptr++) != -1) {
            int s = *ptr++;
            if (!Z[x]) {
                Z[x] = true;
                str[x] = s;
                Q.push(x);
            }
        }
    }

    return true;
}

/**
 * Attract vertices that are in tangles.
 * Current attracting set is <Z>.
 * Current subgame is <Z> + <R>.
 * Write strategy to <str>.
 * Current attracting vertex is <v>.
 * Attracting for player <pl>.
 * (for trace) Attracting to region with priority <pr>.
 */
inline int
RTLSolver::attractTangles(const int pr, const int pl, int v, bitset &R, int *str, bitset &Z)
{
    int added = 0;
    const auto &in_cur = tin[v];
    for (int from : in_cur) {
        if (attractTangle(from, pl, R, str, Z)) {
            added++;
#ifndef NDEBUG
            // maybe report event
            if (trace >= 1) {
                logger << "\033[1;37mattracted \033[1;36m" << tpr[from] << "\033[m-tangle " << from << " to \033[1;36m" << pr << "\033[m" << std::endl;
            }
#endif
        }
    }
    (void)pr; // suppress warning if NDEBUG is defined
    return added;
}

/**
 * Attract vertices that are in tangles.
 * Current attracting set is <Z>.
 * Current subgame is <Z> + <R>.
 * Write strategy to <str>.
 * Attracting for player <pl>.
 * (for trace) Attracting to region with priority <pr>.
 */
inline int
RTLSolver::attractTangles(const int pr, const int pl, bitset &R, int *str, bitset &Z)
{
    int added = 0;
    const unsigned int count = tpr.size();
    for (unsigned int t=0; t<count; t++) {
        if (attractTangle(t, pl, R, str, Z)) {
            added++;
#ifndef NDEBUG
            // maybe report event
            if (trace >= 1) {
                logger << "\033[1;37mattracted \033[1;36m" << tpr[t] << "\033[m-tangle " << t << " to \033[1;36m" << pr << "\033[m" << std::endl;
            }
#endif
        }
    }
    (void)pr; // suppress warning if NDEBUG is defined
    return added;
}

void
RTLSolver::search(bitset &R, int top, int only_player, int depth)
{
    bitset Z(nodecount()); // current region // TODO use the stack thing to only alloc once...

#ifndef NDEBUG
    assert(only_player == -1);
    while (top >= 0 and !R[top]) top--;
    if (top == -1) LOGIC_ERROR; // do not call rec on empty subgame please
#endif

    while (true) {
        // get (next?) top
        while (top >= 0 and !R[top]) top--;
        if (top == -1) return;

        // get top priority and player
        const int pr = priority(top);
        const int pl = pr&1;

        /**
         * Attract from vertices of priority <pr>.
         * Attracted vertices go in Z and are removed from R.
         */

        // We are doing on-the-fly compression (MV trick)

        while (top >= 0 and (priority(top)&1) == (pr&1)) {
            if (R[top] /* and !Z[top]*/) { 
                // if (priority(top) != pr) logger << "DOING IT" << std::endl;
                if (priority(top) != pr) break; // TODO remove this line of course
                heads.push(top); // add to <heads>
                Z[top] = true; // add to <Z> (removed from <R> below)
                str[top] = -1;
                Q.push(top);
            
                while (Q.nonempty()) {
                    const int v = Q.pop();
                    Zvec.push(v);
                    R[v] = false; // remove from R!
                    region[v] = pr; // for debugging mostly
                    attractVertices(pr, pl, v, R, str, Z);
                    attractTangles(pr, pl, v, R, str, Z);
                }
            }
            top--;
        }

        /**
         * If we're not treating this player, report and continue next region...
         */

        if (only_player != -1 and only_player != pl) {
#ifndef NDEBUG
            if (trace >= 2) {
                // report region
                logger << "\033[1;33mregion\033[m ";
                for (int i=0; i<depth; i++) logger << "*";
                logger << "\033[1;36m" << pr << "\033[m";
                for (unsigned int i=0; i<Zvec.size(); i++) {
                    int v = Zvec[i];
                    logger << " \033[1;38;5;160m" << label_vertex(v) << "\033[m";
                    if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                }
                logger << std::endl;
            }
#endif

            // reset for next region
            Z.reset();
            Zvec.clear();
            heads.clear();
            continue;
        }

        /**
         * Identify open heads.
         * Attract for opponent to H and remove from Z.
         */

        unsigned long h_size = 0;
        for (unsigned int i=0; i<heads.size(); i++) {
            int h = heads[i];
            if (!Z[h]) continue; // already removed from Z, ignore
            // check if head is open
            if (owner(h) == pl) {
                if (str[h] != -1) continue; // not open
            } else {
                bool open = false;
                for (auto curedge = outs(h); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (R[to]) { // note: R := R-Z
                        open = true;
                        break;
                    }
                }
                if (!open) continue; // not open
            }
            // head is open, attract for opponent
            H[h] = true; // add to <H> (removed from <Z> below)
            Q.push(h);
            while (Q.nonempty()) {
                const int v = Q.pop();
                h_size++;
                Z[v] = false; // remove from Z!
                attractVertices(pr, 1-pl, v, Z, str, H);
                attractTangles(pr, 1-pl, v, Z, str, H);
            }
        }

        /**
         * If the region is open, then also attract "free" tangles.
         */

        if (h_size != 0) {
            attractTangles(pr, 1-pl, Z, str, H);
            while (Q.nonempty()) {
                const int v = Q.pop();
                h_size++;
                Z[v] = false; // remove from Z!
                attractVertices(pr, 1-pl, v, Z, str, H);
                attractTangles(pr, 1-pl, v, Z, str, H);
            }
        }

#ifndef NDEBUG
        if (trace >= 2) {
            // report region
            logger << "\033[1;33mregion\033[m ";
            for (int i=0; i<depth; i++) logger << "*";
            logger << "\033[1;36m" << pr << "\033[m";
            for (unsigned int i=0; i<Zvec.size(); i++) {
                int v = Zvec[i];
                if (Z[v]) logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                else logger << " \033[1;38;5;160m" << label_vertex(v) << "\033[m";
                if (str[v] != -1) logger << "->" << label_vertex(str[v]);
            }
            logger << std::endl;
        }
#endif

        /**
         * If the region is open, go recursive.
         * Otherwise, compute bottom SCCs as new tangles.
         */

        if (h_size != 0) {
            // reset H
            H.reset();
            // go recursive
            if (Z.any()) {
                heads.clear();
                int _top = Zvec[0];
                Zvec.clear();
                search(Z, _top, only_player, depth+1);
                R &= G;
            }
        } else {
            // extract tangles from each head
            // (because bottom SCC must contain head vertex)

            memset(tarj, 0, sizeof(int[nodecount()]));
            pre = 0;

            while (heads.nonempty()) {
                unsigned int h = heads.pop();
                if (tarj[h] != 0) continue; // already tarjanned
                extractTangles(h, Z, str);
            }

            while (Q.nonempty()) {
                const int v = Q.pop();
                oink->solve(v, pl, str[v]);
                G[v] = false; // remove from Game!
                R[v] = false; // remove from region!
                attractVertices(pr, pl, v, G, str, S);
                attractTangles(pr, pl, v, G, str, S);
            }

            S.reset(); // no need for this anymore.
        }

        // reset for next region
        Z.reset();
        Zvec.clear();
        heads.clear();
    }
}

void
RTLSolver::extractTangles(int i, bitset &R, int *str)
{
    const int pr = priority(i);
    const int pl = pr&1;

    std::vector<int> tangle; // stores the vertices of the tangle

    // start the tarjan search
    Qtar.push(i);
    while (Qtar.nonempty()) {
tarjan_again:
        const unsigned int n = Qtar.back();
        int min; // lowest successor vertex measure
        if (tarj[n] == 0) {
            // first time we see it
            // assign next <pre> to it and put it in <tarres>
            min = (tarj[n] = ++pre);
            tarres.push(n);
        } else {
            min = tarj[n];
        }

        /**
         * Perform the search step of Tarjan.
         */

        if (owner(n) != pl) {
            for (auto curedge = outs(n); *curedge != -1; curedge++) {
                int to = *curedge;
                if (!R[to]) continue; // not in the tangle
                if (tarj[to] == 0) {
                    // not visited, add to search stack and go again
                    Qtar.push(to);
                    goto tarjan_again;
                } else {
                    // visited, update min
                    if (tarj[to] < min) min = tarj[to];
                }
            }
        } else {
            const int s = str[n];
            if (tarj[s] == 0) {
                // not visited, add to search stack and go again
                Qtar.push(s);
                goto tarjan_again;
            } else {
                // visited, update min
                if (tarj[s] < min) min = tarj[s];
            }
        }

        Qtar.pop();

        /**
         * If we're here, then we pushed no new edge.
         * Check if <n> is the root of an SCC.
         */

        if (min < tarj[n]) {
            // Not the root! Update measure and go up.
            tarj[n] = min;
            continue;
        }

        /**
         * Now <n> is the root of an SCC.
         * Every vertex after <n> in <tarres> is also in this SCC.
         * Also set the value of all vertices in the tangle to min
         */

        for (;;) {
            const unsigned int m = tarres.pop();
            tangle.push_back(m);
            tarj[m] = min; // to compute tangleto
            if (m == n) break; // end reached
        }

        /**
         * Now <tangle> contains the current SCC.
         * Every vertex in <tangle> has measure <min>.
         * It is only a real tangle if it contains a cycle, i.e., either 2+ vertices or a self-loop...
         * TODO: if self-loops are removed, then we only need to look at the size!!
         */

        bool is_tangle = (tangle.size() > 1) or
            ((unsigned int)str[n] == n) or
            (str[n] == -1 and game->has_edge(n, n));
        if (!is_tangle) {
            tangle.clear();
            continue;
        }

        /**
         * We have a tangle. Compute the outgoing edges (into <tangleto>) and the next highest region.
         */

#ifndef NDEBUG
        int esc = -1; // escape node (-1 for none, -2 for more)
        int best_esc = -1; // best escape
#endif
        bool bottom_scc = true;
        const auto tangle_end = tangle.end();
        for (auto titer = tangle.begin(); titer != tangle_end;) {
            int v = *titer++;
            if (str[v] != -1) continue; // not losing!
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                int to = *curedge;
                if (disabled[to]) continue; // disabled for solver
                if (tarj[to] == min) continue; // in the tangle
                if (bs_exits[to]) continue; // already added
                bs_exits[to] = true;
                tangleto.push(to);
                /*
                 * we check for bottom SCC's,
                 * because non-bottom SCC's might be duplicate tangles
                 */
                if (tarj[to] != 0) {
                    bottom_scc = false;
                    break;
                }
#ifndef NDEBUG
                // set escape vertex for logging
                if (!bs_exits[v]) {
                    if (esc == -1) esc = v;
                    else if (esc < 0) esc--;
                    else esc = -2;
                    bs_exits[v] = true;
                }
                if (best_esc == -1) best_esc = region[to];
                else if (region[to] < best_esc) best_esc = region[v];
#endif
            }
            if (bottom_scc == false) break;
        }

        // Unmark exits
        bs_exits.reset();

        /**
         * If the tangle is not a bottom SCC, then we skip it.
         * It may have been found already and we don't want duplicates.
         * If it would be a unique tangle, then we'll find it next round!
         */

        if (!bottom_scc) {
            tangle.clear();
            tangleto.clear();
            continue;
        }

        /**
         * If there are no outgoing edges, then we have found a dominion.
         */

        if (tangleto.empty()) {
            // dominion
            if (trace) logger << "\033[1;38;5;201mdominion \033[36m" << pr << "\033[m" << std::endl;
            for (auto titer = tangle.begin(); titer != tangle_end;) {
                int v = *titer++;
                // only if not already solved
                if (S[v] == false) {
                    // actually add to Q and set S[v]
                    Q.push(v);
                    S[v] = true;
                }
            }
            dominions++;
            tangle.clear();
            tangleto.clear();
            continue;
        }

        /**
         * We're not a dominion, we're a tangle.
         */

        /** 
         * First we sort.
         */
        std::sort(tangle.begin(), tangle.end());
        int tsize = tangle.size();

        /**
         * Then we check if it already exists.
         */
        bool already_exists = false;
        for (auto &_tv : tv) {
            for (int i=0; ; i++) {
                if (i == tsize) {
                    if (_tv[2*i] == -1) already_exists=true;
                    break;
                } else if (_tv[2*i] != tangle[i]) {
                    break;
                }
            }
        }
        if (already_exists) {
            if (trace >= 1) {
                logger << "\033[1;38;5;198mduplicate tangle\033[m";
#ifndef NDEBUG
                if (trace >= 2) {
                    for (int &v : tangle) {
                        logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                        if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                    }
                }
#endif
                logger << std::endl;
            }
            tangle.clear();
            tangleto.clear();
            continue;
        }

#ifndef NDEBUG
        if (trace >= 1) {
            logger << "\033[1;38;5;198mnew tangle " << pr << "\033[m (" << tpr.size() << ")";
            if (trace >= 2) {
                for (auto titer = tangle.begin(); titer != tangle_end;) {
                    int v = *titer++;
                    int s = str[v];
                    logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                    if (s != -1) logger << "->" << label_vertex(s);
                }
            }

            logger << " with ";
            if (esc >= 0) logger << "1 escape and ";
            else logger << (-esc) << " escapes and ";
            logger << tangleto.size() << " exit(s)";
            logger << " to region " << best_esc;

            logger << std::endl;
        }
#endif
        // add back links to all normal vertices in our [out]
        const int tidx = tpr.size();
        for (unsigned int x = 0; x < tangleto.size(); x++) tin[tangleto[x]].push_back(tidx);
        // move tangleto into vout
        int* _tout = new int[tangleto.size()+1];
        std::copy(&tangleto[0], &tangleto[tangleto.size()], _tout);
        _tout[tangleto.size()] = -1;
        tout.push_back(_tout);
        // move tangle into vv
        int* _vv = new int[tangle.size()*2+1], c=0;
        for (auto titer = tangle.begin(); titer != tangle_end;) {
            int v = _vv[c++] = *titer++;
            _vv[c++] = str[v];
        }
        _vv[c] = -1;
        tv.push_back(_vv);
        // and set p to pr and current region to BOT
        tpr.push_back(pr);

        tangles++;
        tangleto.clear();
        tangle.clear();
        continue;
    }

    tarres.clear();
}

void
RTLSolver::loop(int player)
{
    bitset R(nodecount());

    while (G.any()) {
        iterations++;

        // store current counts
        const auto D = dominions;
        const auto T = tangles;

        // we do not have any progress yet
        if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "\033[m" << " " << tangles << std::endl;
        R = G;
        search(R, nodecount()-1, player, 0);

        if (D == dominions and T == tangles) {
            // if no dominions or tangles found, return
            return;
        }
    }
}

void
RTLSolver::run()
{
    // get number of nodes and create and initialize inverse array
    max_prio = game->priority(nodecount()-1);

    iterations = 0;
    dominions = 0;
    tangles = 0;

    tin = new std::vector<int>[nodecount()];
    tarj = new int[nodecount()];
    str = new int[nodecount()];
    val = new int[nodecount()]; // for cascader

    // for rec2
    region = new int[nodecount()];
    H.resize(nodecount());
    S.resize(nodecount());
    G = disabled;
    G.flip();
    open_heads.resize(nodecount());

    Q.resize(nodecount());
    Qtar.resize(nodecount());
    Zvec.resize(nodecount());
    heads.resize(nodecount());
    tarres.resize(nodecount());
    tangleto.resize(nodecount());
    bs_exits.resize(nodecount());

    if (onesided) {
        loop(0);
        logger << "after even: " << dominions << " dominions and " << tangles << " tangles and " << iterations << " iterations." << std::endl;

        const int _d = dominions, _t = tangles, _i = iterations;
        iterations = 0;
        dominions = 0;
        tangles = 0;
        // delete all tangles
        for (auto &x : tv) delete[] x;
        for (auto &x : tout) delete[] x;
        tv.clear();
        tout.clear();
        tpr.clear();
        delete[] tin;
        tin = new std::vector<int>[nodecount()];

        loop(1);
        logger << "after odd: " << dominions << " dominions and " << tangles << " tangles and " << iterations << " iterations." << std::endl;

        dominions += _d;
        tangles += _t;
        iterations += _i;
        logger << "found " << dominions << " dominions." << std::endl;
        logger << "solved with " << tangles << " tangles and " << iterations << " iterations." << std::endl;
    } else {
        loop(-1);
        logger << "found " << dominions << " dominions." << std::endl;
        logger << "solved with " << tangles << " tangles and " << iterations << " iterations." << std::endl;
    }

    // check if actually all solved
#ifndef NDEBUG
    for (int i=0; i<nodecount(); i++) {
        if (!disabled[i]) { logger << "incomplete" << std::endl; exit(-1); }
    }
#endif

    // delete[] tangles
    for (auto &x : tv) delete[] x;
    for (auto &x : tout) delete[] x;
    delete[] tin;
    // delete[] tarjan array
    delete[] tarj;
    delete[] str;
}

}
