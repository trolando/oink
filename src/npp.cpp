/*
 * Copyright 2017-2018 Massimo Benerecetti, Fabio Mogavero, Daniele Dell'Erba
 * Copyright 2018 Tom van Dijk
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

#include "npp.hpp"

namespace pg {

NPPSolver::NPPSolver(Oink& oink, Game& game) :
    Solver(oink, game),
    totqueries(0), totpromos(0), maxqueries(0), maxpromos(0), queries(0), promos(0), doms(0),
    maxprio(priority(nodecount() - 1)), strategy(game.getStrategy()), inverse(new int[maxprio + 1]),
    Top(0), End(0), Pivot(0)
{
    // TODO: rewrite to no longer copy the current game.strategy
    uint resprio = maxprio / 20;
    resprio = (resprio >= 500) ? resprio : 500;
    outgame.resize(nodecount());
    winzero.resize(nodecount());
    Phase.reserve(resprio);
    Supgame.reserve(resprio);
    Heads.reserve(resprio);
    Exits.reserve(resprio);
    Entries.reserve(resprio);
    R.resize(nodecount());
    T.resize(nodecount());
    E.resize(nodecount());
}

NPPSolver::~NPPSolver()
{
    delete[] inverse;
    delete Heads[0];
    delete Exits[0];
    delete Entries[0];
    for (uint i = 1; i <= End; ++i)
    {
        delete Supgame[i];
        delete Heads[i];
        delete Exits[i];
        delete Entries[i];
    }
}

_INLINE_ bool NPPSolver::isplayerclosed(uint pos)
{
    for (auto curedge = outs(pos); *curedge != -1; curedge++) {
        int to = *curedge;
        if (R[to])
        {
            strategy[pos] = (!R[pos] || strategy[pos] == -1) ? to : strategy[pos];
            return (true);
        }
    }
    return (false);
}

_INLINE_ bool NPPSolver::isplayerclosedpromo(uint pos)
{
    for (auto curedge = outs(pos); *curedge != -1; curedge++)
    {
        int to = *curedge;
        if (R[to])
        {
            strategy[pos] = to;
            return (true);
        }
    }
    return (false);
}

_INLINE_ bool NPPSolver::isopponentclosedongame(uint pos)
{
    for (auto curedge = outs(pos); *curedge != -1; curedge++)
    {
        int to = *curedge;
        if (outgame[to] || R[to])
        {
            continue;
        }
        else
        {
            return (false);
        }
    }
    strategy[pos] = -1;
    return (true);
}

_INLINE_ bool NPPSolver::isopponentclosedonsubgame(uint pos)
{
    auto & supgame = *(Supgame[Top]);
    for (auto curedge = outs(pos); *curedge != -1; curedge++)
    {
        int to = *curedge;
        if (outgame[to] || R[to])
        {
            continue;
        }
        else if (supgame[to])
        {
            E.push(to);
        }
        else
        {
            E.clear();
            return (false);
        }
    }
    auto & exits = *(Exits[Pivot]);
    while (E.nonempty())
    {
        exits[E.pop()] = true;
    }
    strategy[pos] = -1;
    return (true);
}

_INLINE_ bool NPPSolver::isclosedongame(uint pos)
{
    if ((uint) owner(pos) == alpha)
    {
        return (isplayerclosed(pos));
    }
    else
    {
        return (isopponentclosedongame(pos));
    }
}

_INLINE_ bool NPPSolver::isclosedonsubgame(uint pos)
{
    if ((uint) owner(pos) == alpha)
    {
        return (isplayerclosed(pos));
    }
    else
    {
        return (isopponentclosedonsubgame(pos));
    }
}

_INLINE_ bool NPPSolver::isclosedonsubgamepromo(uint pos)
{
    if ((uint) owner(pos) == alpha)
    {
        return (isplayerclosedpromo(pos));
    }
    else
    {
        return (isopponentclosedonsubgame(pos));
    }
}

_INLINE_ void NPPSolver::pushinqueueongame(uint pos)
{
    for (auto curedge = ins(pos); *curedge != -1; curedge++)
    {
        int from = *curedge;
        if (outgame[from] || R[from])
        {
            continue;
        }
        else if ((uint) owner(from) == alpha)
        {
            R[from] = true;
            strategy[from] = pos;
            T.push(from);
        }
        else
        {
            if (isopponentclosedongame(from))
            {
                R[from] = true;
                T.push(from);
            }
        }
    }
}

bool NPPSolver::atrongame()
{
    auto & entries = *(Entries[Top]);
    auto lend = entries.end();
    for (auto liter = entries.begin(); liter != lend; ++liter)
    {
        auto & rdeque = *liter;
        auto dend = rdeque.end();
        for (auto diter = rdeque.begin(); diter != dend; ++diter)
        {
            pos = *diter;
            if (!R[pos] && isclosedongame(pos))
            {
                R[pos] = true;
                T.push(pos);
            }
        }
    }
    if (T.nonempty())
    {
        do
        {
            pushinqueueongame(T.pop());
        }
        while (T.nonempty());
        return (true);
    }
    else
    {
        return (false);
    }
}

_INLINE_ void NPPSolver::pushinqueueonsubgamedw(uint pos)
{
    auto & supgame = *(Supgame[Top]);
    auto & entries = *(Entries[Top]->begin());
    for (auto curedge = ins(pos); *curedge != -1; curedge++)
    {
        int from = *curedge;
        if (outgame[from])
        {
            continue;
        }
        else if (supgame[from])
        {
            entries.push_back(from);
            continue;
        }
        else if ((uint) owner(from) == alpha)
        {
            if (!R[from])
            {
                R[from] = true;
                strategy[from] = pos;
                T.push(from);
            }
            else if (O[from])
            {
                O[from] = false;
                strategy[from] = pos;
            }
        }
        else if (!R[from])
        {
            if (isopponentclosedonsubgame(from))
            {
                R[from] = true;
                T.push(from);
            }
        }
        else if (O[from])
        {
            if (isopponentclosedonsubgame(from))
            {
                O[from] = false;
            }
        }
    }
}

void NPPSolver::atronsubgamedw()
{
    auto & heads = *(Heads[Top]);
    auto qend = heads.end();
    for (auto qiter = heads.begin(); qiter != qend; ++qiter)
    {
        pushinqueueonsubgamedw(*qiter);
    }
    while (T.nonempty())
    {
        pushinqueueonsubgamedw(T.pop());
    }
}

_INLINE_ void NPPSolver::pushinqueueonsubgameup(uint pos)
{
    auto & supgame = *(Supgame[Top]);
    auto & entries = *(Entries[Pivot]->begin());
    for (auto curedge = ins(pos); *curedge != -1; curedge++)
    {
        int from = *curedge;
        if (outgame[from] || R[from])
        {
            continue;
        }
        else if (supgame[from])
        {
            entries.push_back(from);
            continue;
        }
        else if ((uint) owner(from) == alpha)
        {
            R[from] = true;
            strategy[from] = pos;
            T.push(from);
        }
        else
        {
            if (isopponentclosedonsubgame(from))
            {
                R[from] = true;
                T.push(from);
            }
        }
    }
}

bool NPPSolver::atronsubgameup()
{
    auto & supgame = *(Supgame[Top]);
    auto & entries = *(Entries[Pivot]);
    auto lend = entries.end();
    for (auto liter = entries.begin(); liter != lend; ++liter)
    {
        auto & rdeque = *liter;
        auto dend = rdeque.end();
        for (auto diter = rdeque.begin(); diter != dend; ++diter)
        {
            pos = *diter;
            if (!supgame[pos] && !R[pos] && isclosedonsubgame(pos))
            {
                R[pos] = true;
                T.push(pos);
            }
        }
    }
    if (T.nonempty())
    {
        do
        {
            pushinqueueonsubgameup(T.pop());
        }
        while (T.nonempty());
        return (true);
    }
    else
    {
        return (false);
    }
}

_INLINE_ void NPPSolver::goup()
{
    D = *(Supgame[Top]) - *(Supgame[Top - 1]);
    p = priority(Heads[--Top]->front());
    beta = p & 1;
}

_INLINE_ void NPPSolver::newstackslot()
{
    Phase.push_back(true);
    Supgame.push_back(new bitset(R | *(Supgame[Top])));
    Pivot = End = ++Top;
    Heads.push_back(new uideque());
    Exits.push_back(new bitset(nodecount()));
    Entries.push_back(new uidlist());
    Entries[Top]->push_front(uideque());
}

_INLINE_ void NPPSolver::clearstackslot()
{
    Pivot = ++Top;
    Phase[Top] = true;
    *(Supgame[Top]) = *(Supgame[Top - 1]);
    *(Supgame[Top]) |= R;
    Heads[Top]->clear();
    Exits[Top]->reset();
    Entries[Top]->clear();
    Entries[Top]->push_front(uideque());
}

_INLINE_ void NPPSolver::nextpriopos()
{
    auto & supgame = *(Supgame[Top]);
    int pstar;
    for (pstar = p - 1; inverse[pstar] == -1; --pstar)
    {
    }
    for (p = pstar, pos = inverse[p]; supgame[pos]; p = priority(--pos))
    {
    }
    alpha = p & 1;
    R.reset();
}

_INLINE_ void NPPSolver::godwdw()
{
    Phase[Top] = false;
    if (Top >= End) // A new stack slot needs to be allocated
    {
        newstackslot();
    }
    else // It can reuse some previous stack slot
    {
        clearstackslot();
    }
    nextpriopos();
}

_INLINE_ void NPPSolver::godwup()
{
    Phase[Top] = false;
    clearstackslot();
    nextpriopos();
}

_INLINE_ void NPPSolver::search()
{

    /* vv Resetting local statistics vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
    queries = 0;
    promos = 0;
    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

    while (true)
    {
        if (Phase[Top])
        {

            /* vv Update of the statistic on the number of queries vvvvvvvvvvvvvvvv */
            ++queries;
            /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

            /* vv Search of region heads vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
            auto & supgame = *(Supgame[Top]);
            auto & heads = *(Heads[Top]);
            for (; pos >= 0 && (uint) priority(pos) == p; --pos) // Collect heads of the region
            {
                if (!supgame[pos])
                {
                    R[pos] = true;
                    heads.push_front(pos);
                }
            }
            O = R; // Mark the heads of the new region as open heads
            /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

            /* vv Region construction via attractor vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
            atronsubgamedw();
            /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

            /* vv Determine the successor state vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
            if (O.none()) // No heads of the region are open
            {
                if (Exits[Top]->none())
                {
                    break; // Region closed in the whole game (i.e. it is dominion)
                }
                else
                {
                    goup(); // Region closed in the subgame, look upward for the target of the promotion
                }
            }
            else
            {
                godwdw(); // Region open, descend further
            }
            /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

        }
        else
        {
            auto & pivotexits = *(Exits[Pivot]);
            if (alpha != beta || (D & pivotexits).none())
            {
                goup(); // Target D is of the opponent or R cannot enter D (i.e., BEP not found yet), go further up
            }
            else // Target of the promotion found (BEP)
            {

                /* vv Update of the statistics on the number of queries and promotions*/
                ++queries;
                ++promos;
                /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

                /* vv Region merging and maximization via attractor vvvvvvvvvvvvvvvvv */
                R |= D;
                atronsubgameup();
                /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

                /* vv Closure check vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
                bool closed = true;
                auto & heads = *(Heads[Top]);
                while (!heads.empty())
                {
                    pos = heads.front();
                    if (isclosedonsubgamepromo(pos))
                    {
                        heads.pop_front();
                    }
                    else
                    {
                        closed = false;
                        break;
                    }
                }
                /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

                /* vv Update of promotion exists vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
                pivotexits -= D;
                *(Exits[Top]) |= pivotexits;
                /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

                /* vv Update of potential region entries vvvvvvvvvvvvvvvvvvvvvvvvvvvv */
                auto & topentries = *(Entries[Top]);
                auto & pivotentries = *(Entries[Pivot]);
                topentries.splice(topentries.end(), pivotentries);
                /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

                /* vv Determine the successor state vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
                Pivot = Top;
                if (closed)
                {
                    if (Exits[Top]->none())
                    {
                        break; // Region closed in the game (i.e. it is dominion)
                    }
                    else
                    {
                        goup(); // Region closed in the subgame, look upward for the target of the promotion
                    }
                }
                else
                {
                    godwup(); // Region open, start new descent
                }
                /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

            }
        }
    }

    /* vv Update of the statistics on the number of queries and promotions vvvv */
    totqueries += queries;
    totpromos += promos;
    if (queries > maxqueries)
    {
        maxqueries = queries;
    }
    if (promos > maxpromos)
    {
        maxpromos = promos;
    }
    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

    /* vv Dominion extension via attractor vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
    atrongame();
    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

}

void NPPSolver::run()
{

    /* vv Stack initialization vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
    Phase.push_back(true);
    Supgame.push_back(&outgame);
    Heads.push_back(new uideque());
    Exits.push_back(new bitset(nodecount()));
    Entries.push_back(new uidlist());
    Entries[0]->push_front(uideque());
    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

    /* vv Initialization of pos and maxprio vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
    for (pos = nodecount() - 1;; --pos)
    {
        if (disabled[pos])
        {
            outgame[pos] = true;
        }
        else
        {
            maxprio = priority(pos);
            break;
        }
    }
    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

    /* vv Initialization of inverse, outgame, and strategy vvvvvvvvvvvvvvvvvvvv */
    for (uint prt = 0; prt <= maxprio; ++prt)
    {
        inverse[prt] = -1;
    }
    for (int sop = 0; sop <= pos; sop++)
    {
        if (disabled[sop])
        {
            outgame[sop] = true;
        }
        else
        {
            inverse[priority(sop)] = sop;
            strategy[sop] = -1;
        }
    }
    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

    /* vv Main solution cycle vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
    while (pos >= 0)
    {

        /* vv Update of the statistic on the number of dominions vvvvvvvvvvvvvvvv */
        doms++;
        /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

        /* vv Setting of the initial priority and parity vvvvvvvvvvvvvvvvvvvvvvvv */
        p = maxprio;
        alpha = p & 1;
        /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

        /* vv Call to the search routine and update of the zero winning region vv */
        search();
        outgame |= R;
        if (alpha == 0)
        {
            winzero |= R;
        }
        /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

        /* vv Search for the new pos and maxprio vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
        for (pos = inverse[maxprio]; pos >= 0 && outgame[pos]; --pos)
        {
        }
        if (pos >= 0) maxprio = priority(pos);
        /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

        /* vv Stack and current region cleaning vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
        Top = Pivot = 0;
        Phase[0] = true;
        Heads[0]->clear();
        Exits[0]->reset();
        Entries[0]->clear();
        Entries[0]->push_front(uideque());
        R.reset();
        // /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

    }
    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

    /* vv Setting of the final solution vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
    for (pos = 0; pos < nodecount(); ++pos)
    {
        if (disabled[pos]) continue;
        Solver::solve(pos, !winzero[pos], strategy[pos]);
    }
    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */

    /* vv Printing of statistics vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv */
    logger << "solved with " << totqueries << " total queries and " << totpromos << " total promotions;" << std::endl;
    logger << "            " << maxqueries << " max queries and " << maxpromos << " max promotions;" << std::endl;
    logger << "            " << doms << " dominions." << std::endl;
    /* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
}

}
