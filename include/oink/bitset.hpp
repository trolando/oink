/**
 * Copyright 2020 Tom van Dijk
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

#ifndef BITSET_HPP
#define BITSET_HPP

#include <oink/libpopcnt.h>

namespace pg
{

static inline int bsr(uint64_t x)
{
    return __builtin_clzll(x) ^ 63;
}

class bitset
{
public:
    class reference
    {
        friend class bitset;

        reference(uint64_t &b, unsigned int pos) : _block(b), _mask(uint64_t(1)<<pos) { }
        reference(const reference& other) = delete;

        void operator&(); // left undefined

    public:
        operator bool() const { return (_block & _mask) != 0; }
        bool operator~() const { return (_block & _mask) == 0; }

        reference& flip() { do_flip(); return *this; }

        reference& operator=(bool x)               { do_assign(x);   return *this; } // for b[i] = x
        reference& operator=(const reference& rhs) { do_assign(rhs); return *this; } // for b[i] = b[j]

        reference& operator|=(bool x) { if  (x) do_set();   return *this; }
        reference& operator&=(bool x) { if (!x) do_reset(); return *this; }
        reference& operator^=(bool x) { if  (x) do_flip();  return *this; }
        reference& operator-=(bool x) { if  (x) do_reset(); return *this; }

     private:
        uint64_t &_block;
        const uint64_t _mask;

        void do_set() { _block |= _mask; }
        void do_reset() { _block &= ~_mask; }
        void do_flip() { _block ^= _mask; }
        void do_assign(bool x) { x ? do_set() : do_reset(); }
    };

    bitset()
    {
        _size = 0;
        _bitssize = 0;
        _allocsize = 0;
        _bits = NULL;
    }

    bitset(size_t newsize)
    {
        _size = newsize;
        _bitssize = (_size+63)/64;
        _allocsize = _bitssize * 8;
        _bits = new uint64_t[_bitssize];
        std::fill(_bits, _bits+_bitssize, '\0');
    }

    bitset(const bitset &other)
    {
        _size = other._size;
        _bitssize = other._bitssize;
        _allocsize = _bitssize * 8;
        _bits = new uint64_t[_bitssize];
        std::copy(other._bits, other._bits+_bitssize, _bits);
    }

    ~bitset()
    {
        if (_bits != NULL) delete[] _bits;
    }

    /**
     * After resizing the *new* bits _can_ be undefined.
     */
    void resize(size_t newsize)
    {
        if (_allocsize == 0) {
            bitset b(newsize);
            swap(b);
        } else if (newsize <= _allocsize*8) {
            // fits in allocated array already
            _size = newsize;
            _bitssize = (newsize+63)/64;
            zero_unused_bits(); // only zeroes the last used block...
        } else {
            bitset b(newsize);
            std::copy(_bits, _bits+_bitssize, b._bits);
            swap(b);
        }
    }

private:
    inline size_t num_blocks(void) const { return _bitssize; }
    inline size_t block_index(size_t pos) const { return pos / 64; }
    inline size_t bit_index(size_t pos) const { return pos % 64; }
    inline uint64_t bit_mask(size_t pos) const { return uint64_t(1) << bit_index(pos); }
    inline size_t count_extra_bits(void) const { return _size % 64; }

    inline void zero_unused_bits()
    {
        size_t extra = count_extra_bits();
        if (extra != 0) _bits[num_blocks()-1] &= ((uint64_t(1) << extra) - 1);
    }

public:
    __attribute__((always_inline)) bitset& reset(void)
    {
        std::fill(_bits, _bits+num_blocks(), '\0');
        return *this;
    }

    __attribute__((always_inline)) bitset& set(void)
    {
        std::fill(_bits, _bits+num_blocks(), static_cast<uint64_t>(~0));
        zero_unused_bits();
        return *this;
    }

    inline bitset& flip(void)
    {
        for (size_t i=0; i<num_blocks(); i++) _bits[i] = ~_bits[i];
        zero_unused_bits();
        return *this;
    }

    inline std::size_t size(void) const
    {
        return _size;
    }

    std::size_t count(void) const
    {
        return popcnt(_bits, num_blocks()*8);
    }

    inline bool any(void) const
    {
        uint64_t *p = _bits;
        std::size_t len = num_blocks();
        while (len-- != 0) if (*p++) return true;
        return false;
    }

    inline bool none() const
    {
        return !any();
    }

    inline bool all() const
    {
        if (empty()) return true;
        size_t extra = count_extra_bits();
        if (extra == 0) {
            uint64_t *p = _bits;
            std::size_t len = num_blocks();
            while (len-- != 0) if (*p++ != static_cast<uint64_t>(~0)) return false;
        } else {
            uint64_t *p = _bits;
            std::size_t len = num_blocks()-1;
            while (len-- != 0) if (*p++ != static_cast<uint64_t>(~0)) return false;
            const uint64_t last_mask = (uint64_t(1)<<extra)-1;
            if (*p != last_mask) return false;
        }
        return true;
    }

    inline bool empty() const
    {
        return _size == 0;
    }

    bitset operator~() const
    {
        bitset b(*this);
        b.flip();
        return b;
    }

    inline void reset(size_t pos)
    {
        _bits[block_index(pos)] &= ~bit_mask(pos);
    }

    inline void set(size_t pos)
    {
        _bits[block_index(pos)] |= bit_mask(pos);
    }

    inline bool test(size_t pos) const
    {
        return (_bits[block_index(pos)] & bit_mask(pos)) != 0;
    }

    reference operator[](size_t pos)
    {
        return reference(_bits[block_index(pos)], bit_index(pos));
    }

    inline bool operator[](size_t pos) const
    {
        return test(pos);
    }

    bitset& operator=(const bitset &src)
    {
        bitset b(src);
        swap(b);
        return *this;
    }

    bitset& operator-=(const bitset& other)
    {
        uint64_t *p = _bits;
        const uint64_t *q = other._bits;
        std::size_t len = num_blocks();
        while (len-- != 0) (*p++) &= ~(*q++);
        return *this;
    }

    bitset& operator&=(const bitset& other)
    {
        uint64_t *p = _bits;
        const uint64_t *q = other._bits;
        std::size_t len = num_blocks();
        while (len-- != 0) (*p++) &= (*q++);
        return *this;
    }

    bitset& operator|=(const bitset &other)
    {
        uint64_t *p = _bits;
        const uint64_t *q = other._bits;
        std::size_t len = num_blocks();
        while (len-- != 0) (*p++) |= (*q++);
        return *this;
    }

    bitset& operator^=(const bitset &other)
    {
        uint64_t *p = _bits;
        const uint64_t *q = other._bits;
        std::size_t len = num_blocks();
        while (len-- != 0) (*p++) ^= (*q++);
        return *this;
    }

    bool operator==(const bitset &other) const
    {
        const uint64_t *p = _bits;
        const uint64_t *q = other._bits;
        std::size_t len = num_blocks();
        while (len-- != 0) if ((*p++) != (*q++)) return false;
        return true;
    }

    bool operator!=(const bitset &other) const
    {
        return !(*this == other);
    }

    inline void swap(bitset &other)
    {
        std::swap(_size, other._size);
        std::swap(_bitssize, other._bitssize);
        std::swap(_allocsize, other._allocsize);
        std::swap(_bits, other._bits);
    }

    bool intersects(const bitset& other) const
    {
        const uint64_t *p = _bits;
        const uint64_t *q = other._bits;
        std::size_t len = num_blocks();
        while (len-- != 0) if ((*p++) & (*q++)) return true;
        return false;
    }

    size_t find_first() const
    {
        size_t i = 0;
        while (i < num_blocks() and _bits[i] == 0) i++;
        if (i == num_blocks()) return npos;
        else return i*64 + __builtin_ffsll(_bits[i]) - 1;
    }

    size_t find_last() const
    {
        if (num_blocks() == 0) return npos;

        size_t i = num_blocks()-1;
        for (;;) {
            if (_bits[i] != 0) return i*64 + bsr(_bits[i]);
            if (i == 0) return npos;
            i--;
        }
    }

    size_t find_next(size_t pos) const
    {
        if (pos == npos or (pos+1) >= _size) return npos;
        pos++;
        size_t i = block_index(pos);
        uint64_t m = _bits[i] & (~uint64_t(0) << bit_index(pos));
        if (m) {
            return i*64 + __builtin_ffsll(m) - 1;
        } else {
            i += 1;
            while (i < num_blocks() and _bits[i] == 0) i++;
            if (i == num_blocks()) return npos;
            else return i*64 + __builtin_ffsll(_bits[i]) - 1;
        }
    }

    size_t find_prev(size_t pos) const
    {
        if (pos == 0 or pos == npos) return npos;
        size_t i = block_index(pos);
        uint64_t m = _bits[i] & ~((~uint64_t(0)) << (int)bit_index(pos));
        if (m) {
            return i*64 + bsr(m);
        } else {
            if (i == 0) return npos;
            i -= 1;
            for (;;) {
                if (_bits[i] != 0) return i*64 + bsr(_bits[i]);
                if (i == 0) return npos;
                i--;
            }
        }
    }

    static const size_t npos = static_cast<size_t>(-1);

protected:
    uint64_t *_bits;
    size_t _size, _bitssize;
    size_t _allocsize;
};

inline bitset operator^(const bitset& x, const bitset& y)
{
    bitset b(x);
    return b ^= y;
}

inline bitset operator-(const bitset& x, const bitset& y)
{
    bitset b(x);
    return b -= y;
}

inline bitset operator|(const bitset& x, const bitset &y)
{
    bitset b(x);
    return b |= y;
}

inline bitset operator&(const bitset& x, const bitset &y)
{
    bitset b(x);
    return b &= y;
}

inline void swap(bitset &left, bitset &right)
{
    left.swap(right);
}

}

namespace std
{
    inline void swap(pg::bitset &left, pg::bitset &right)
    {
        left.swap(right);
    }
}

#endif
