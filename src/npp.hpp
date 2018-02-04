/*
 * Copyright 2017-2018 Massimo Benerecetti, Fabio Mogavero, Daniele Dell'Erba
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

#ifndef NPP_HPP
#define NPP_HPP

#include <list>
#include <deque>
#include <vector>
#include <boost/dynamic_bitset.hpp>

#include "oink.hpp"
#include "solver.hpp"

#define _INLINE_ __attribute__((always_inline))

namespace pg
{

class uintqueue
{
public:

    /******************************************************************************/
    /* Typedef                                                                    */
    /******************************************************************************/

    typedef unsigned int uint;

    /******************************************************************************/

    /******************************************************************************/
    /* Constructors and destructor                                                */
    /******************************************************************************/

    uintqueue()
    {
        pointer = 0;
        queue = NULL;
    }

    uintqueue(uint capacity)
    {
        pointer = 0;
        queue = new uint[capacity];
    }

    ~uintqueue()
    {
        delete[] queue;
    }

    /******************************************************************************/

    /******************************************************************************/
    /* Main methods                                                               */
    /******************************************************************************/

    _INLINE_ bool nonempty()
    {
        return ((bool) pointer);
    }

    _INLINE_ bool empty()
    {
        return pointer == 0;
    }

    _INLINE_ uint pop()
    {
        return (queue[--pointer]);
    }

    _INLINE_ void push(uint element)
    {
        queue[pointer++] = element;
    }

    _INLINE_ uint& back()
    {
        return (queue[pointer-1]);
    }

    _INLINE_ void clear()
    {
        pointer = 0;
    }

    _INLINE_ uint& operator[] (uint idx)
    {
        return (queue[idx]);
    }

    _INLINE_ uint size()
    {
        return (pointer);
    }

    _INLINE_ void resize(uint capacity)
    {
        pointer = 0;
        if (queue != NULL)
        {
            delete[] queue;
        }
        queue = new uint[capacity];
    }

    _INLINE_ void swap(uintqueue & other)
    {
        uint * tmp = other.queue;
        other.queue = queue;
        queue = tmp;
    }

    /******************************************************************************/

protected:

    /******************************************************************************/
    /* Queue and pointer                                                          */
    /******************************************************************************/

    uint * queue;
    int pointer;

    /******************************************************************************/

};

class NPPSolver : public Solver
{
public:

    /******************************************************************************/
    /* Constructor and destructor                                                 */
    /******************************************************************************/

    NPPSolver(Oink * oink, Game * game, std::ostream & lgr = std::cout);

    virtual ~NPPSolver();

    /******************************************************************************/

    /******************************************************************************/
    /* Main method                                                                */
    /******************************************************************************/

    // Main solver function
    virtual void run();

    /******************************************************************************/

protected:

    /******************************************************************************/
    /* Typedefs                                                                   */
    /******************************************************************************/

    typedef unsigned int uint;

    typedef std::deque<uint> uideque;

    typedef std::list<uideque> uidlist;

    typedef boost::dynamic_bitset<unsigned long long> bitset;

    typedef std::vector<bool> bvector;

    typedef std::vector<uideque *> uidvector;

    typedef std::vector<bitset *> bsvector;

    typedef std::vector<uidlist *> uidlvector;

    /******************************************************************************/

    /******************************************************************************/
    /* Statistic fields                                                           */
    /******************************************************************************/

    uint totqueries;  // Total number of queries needed for the game solution
    uint totpromos;   // Total number of promotions needed for the game solution

    uint maxqueries;  // Maximal number of queries needed for finding a dominion
    uint maxpromos;   // Maximal number of promotions needed for finding a dominion

    uint queries;     // Current number of queries
    uint promos;      // Current number of promotions

    uint doms;        // Number of dominions found

    /******************************************************************************/

    /******************************************************************************/
    /* Game fields                                                                */
    /******************************************************************************/

    uint maxprio;   // Maximal priority

    int * strategy; // Strategies of the players

    int * inverse;  // Maximal positions associated with the priorities

    bitset outgame; // Positions whose winner has been already determined

    bitset winzero; // Positions already won by player zero

    /******************************************************************************/

    /******************************************************************************/
    /* Stack fields                                                               */
    /******************************************************************************/

    uint Top;           // Index of the current element within the vectors
    uint End;           // Index of the last element within the vectors
    uint Pivot;         // Index of the element associated with the region variables

    bvector Phase;      // Stack of phases

    bsvector Supgame;   // Stack of supgames

    uidvector Heads;    // Stack of region heads

    bsvector Exits;     // Stack of region promotion exits

    uidlvector Entries; // Stack of region potential entriers

    /******************************************************************************/

    /******************************************************************************/
    /* Current region fields                                                      */
    /******************************************************************************/

    bitset R;   // Positions in the current region

    uint alpha; // Player of the current region R

    /******************************************************************************/

    /******************************************************************************/
    /* Previous region fields                                                     */
    /******************************************************************************/

    bitset D;   // Positions in the previous region upward along the stack

    uint beta;  // Previous region player (player of region D)

    /******************************************************************************/

    /******************************************************************************/
    /* Accessory region fields                                                    */
    /******************************************************************************/

    int pos;      // Working position

    uint p;       // Working priority

    bitset O;     // Bitset of open heads of the current region R

    uintqueue T;  // Tail queue for positions waiting for attraction

    uintqueue E;  // Exits queue for positions waiting for inclusion as exits

    /******************************************************************************/

    /******************************************************************************/
    /* Accessory methods                                                          */
    /******************************************************************************/

    // Attractor functions
    virtual bool atrongame();
    virtual void atronsubgamedw();
    virtual bool atronsubgameup();

    /******************************************************************************/

    /******************************************************************************/
    /* Solver method                                                              */
    /******************************************************************************/

    // Search function
    inline virtual void search();

    /******************************************************************************/

private:

    /******************************************************************************/
    /* Accessory methods                                                           */
    /******************************************************************************/

    // Go up/down functions
    inline void goup();
    inline void godwdw();
    inline void godwup();

    // Top stack handler functions
    inline void newstackslot();
    inline void clearstackslot();

    // Next priority and position function
    inline void nextpriopos();

    // Push in queue functions
    inline void pushinqueueongame(uint pos);
    inline void pushinqueueonsubgamedw(uint pos);
    inline void pushinqueueonsubgameup(uint pos);

    // Is closed functions
    inline bool isclosedongame(uint pos);
    inline bool isclosedonsubgame(uint pos);
    inline bool isclosedonsubgamepromo(uint pos);
    inline bool isplayerclosed(uint pos);
    inline bool isplayerclosedpromo(uint pos);
    inline bool isopponentclosedongame(uint pos);
    inline bool isopponentclosedonsubgame(uint pos);

    /******************************************************************************/
};

}

#endif
