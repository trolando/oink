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

#ifndef UINTQUEUE_HPP
#define UINTQUEUE_HPP

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

    __attribute__((always_inline)) bool nonempty()
    {
        return ((bool) pointer);
    }

    __attribute__((always_inline)) bool empty()
    {
        return pointer == 0;
    }

    __attribute__((always_inline)) uint pop()
    {
        return (queue[--pointer]);
    }

    __attribute__((always_inline)) void pop2()
    {
        pointer -= 2;
    }

    __attribute__((always_inline)) void push(uint element)
    {
        queue[pointer++] = element;
    }

    __attribute__((always_inline)) uint& back()
    {
        return (queue[pointer-1]);
    }

    __attribute__((always_inline)) uint& back2()
    {
        return (queue[pointer-2]);
    }

    __attribute__((always_inline)) void clear()
    {
        pointer = 0;
    }

    __attribute__((always_inline)) uint& operator[] (uint idx)
    {
        return (queue[idx]);
    }

    __attribute__((always_inline)) uint size()
    {
        return (pointer);
    }

    __attribute__((always_inline)) void resize(uint capacity)
    {
        pointer = 0;
        if (queue != NULL)
        {
            delete[] queue;
        }
        queue = new uint[capacity];
    }

    __attribute__((always_inline)) void swap(uintqueue & other)
    {
        std::swap(queue, other.queue);
        std::swap(pointer, other.pointer);
    }

    __attribute__((always_inline)) void swap_elements(uint idx1, uint idx2)
    {
        std::swap(queue[idx1], queue[idx2]);
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

}

#endif
