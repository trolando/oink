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

}

#endif
