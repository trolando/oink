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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <oink/intrinsics.hpp>
#include <oink/libpopcnt.h>

namespace pg
{

static inline int bsr(uint64_t x)
{
    return intrinsics::countl_zero64(x) ^ 63;
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

    bitset() : _size(0) { }

    bitset(size_t newsize) : bits_((newsize+63)/64, 0), _size(newsize) { }

    bitset(const bitset &other) = default;

    ~bitset() = default;

    /**
     * After resizing, any new bits are zero.
     */
    void resize(size_t newsize)
    {
        _size = newsize;
        bits_.resize((newsize+63)/64);
        zero_unused_bits();
    }

private:
    inline size_t num_blocks(void) const { return bits_.size(); }
    inline size_t block_index(size_t pos) const { return pos / 64; }
    inline size_t bit_index(size_t pos) const { return pos % 64; }
    inline uint64_t bit_mask(size_t pos) const { return uint64_t(1) << bit_index(pos); }
    inline size_t count_extra_bits(void) const { return _size % 64; }

    inline void zero_unused_bits()
    {
        size_t extra = count_extra_bits();
        if (extra != 0) bits_[num_blocks()-1] &= ((uint64_t(1) << extra) - 1);
    }

public:
    /** Low-level access to the underlying 64-bit blocks. */
    [[nodiscard]] uint64_t* data() noexcept { return bits_.data(); }
    [[nodiscard]] const uint64_t* data() const noexcept { return bits_.data(); }
    [[nodiscard]] std::size_t block_count() const noexcept { return bits_.size(); }

    __attribute__((always_inline)) bitset& reset(void)
    {
        std::fill(bits_.begin(), bits_.end(), uint64_t(0));
        return *this;
    }

    __attribute__((always_inline)) bitset& set(void)
    {
        std::fill(bits_.begin(), bits_.end(), static_cast<uint64_t>(~0));
        zero_unused_bits();
        return *this;
    }

    inline bitset& flip(void)
    {
        for (size_t i=0; i<num_blocks(); i++) bits_[i] = ~bits_[i];
        zero_unused_bits();
        return *this;
    }

    inline std::size_t size(void) const
    {
        return _size;
    }

    std::size_t count(void) const
    {
        return popcnt(bits_.data(), num_blocks()*8);
    }

    inline bool any(void) const
    {
        const uint64_t *p = bits_.data();
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
            const uint64_t *p = bits_.data();
            std::size_t len = num_blocks();
            while (len-- != 0) if (*p++ != static_cast<uint64_t>(~0)) return false;
        } else {
            const uint64_t *p = bits_.data();
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
        bits_[block_index(pos)] &= ~bit_mask(pos);
    }

    inline void set(size_t pos)
    {
        bits_[block_index(pos)] |= bit_mask(pos);
    }

    inline bool test(size_t pos) const
    {
        return (bits_[block_index(pos)] & bit_mask(pos)) != 0;
    }

    reference operator[](size_t pos)
    {
        return reference(bits_[block_index(pos)], bit_index(pos));
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
        uint64_t *p = bits_.data();
        const uint64_t *q = other.bits_.data();
        std::size_t len = num_blocks();
        while (len-- != 0) (*p++) &= ~(*q++);
        return *this;
    }

    bitset& operator&=(const bitset& other)
    {
        uint64_t *p = bits_.data();
        const uint64_t *q = other.bits_.data();
        std::size_t len = num_blocks();
        while (len-- != 0) (*p++) &= (*q++);
        return *this;
    }

    bitset& operator|=(const bitset &other)
    {
        uint64_t *p = bits_.data();
        const uint64_t *q = other.bits_.data();
        std::size_t len = num_blocks();
        while (len-- != 0) (*p++) |= (*q++);
        return *this;
    }

    bitset& operator^=(const bitset &other)
    {
        uint64_t *p = bits_.data();
        const uint64_t *q = other.bits_.data();
        std::size_t len = num_blocks();
        while (len-- != 0) (*p++) ^= (*q++);
        return *this;
    }

    bool operator==(const bitset &other) const
    {
        const uint64_t *p = bits_.data();
        const uint64_t *q = other.bits_.data();
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
        bits_.swap(other.bits_);
        std::swap(_size, other._size);
    }

    bool intersects(const bitset& other) const
    {
        const uint64_t *p = bits_.data();
        const uint64_t *q = other.bits_.data();
        std::size_t len = num_blocks();
        while (len-- != 0) if ((*p++) & (*q++)) return true;
        return false;
    }

    size_t find_first() const
    {
        size_t i = 0;
        while (i < num_blocks() and bits_[i] == 0) i++;
        if (i == num_blocks()) return npos;
        else return i*64 + intrinsics::countr_zero64(bits_[i]);
    }

    size_t find_last() const
    {
        if (num_blocks() == 0) return npos;

        size_t i = num_blocks()-1;
        for (;;) {
            if (bits_[i] != 0) return i*64 + bsr(bits_[i]);
            if (i == 0) return npos;
            i--;
        }
    }

    size_t find_next(size_t pos) const
    {
        if (pos == npos or (pos+1) >= _size) return npos;
        pos++;
        size_t i = block_index(pos);
        uint64_t m = bits_[i] & (~uint64_t(0) << bit_index(pos));
        if (m) {
            return i*64 + intrinsics::countr_zero64(m);
        } else {
            i += 1;
            while (i < num_blocks() and bits_[i] == 0) i++;
            if (i == num_blocks()) return npos;
            else return i*64 + intrinsics::countr_zero64(bits_[i]);
        }
    }

    size_t find_prev(size_t pos) const
    {
        if (pos == 0 or pos == npos) return npos;
        size_t i = block_index(pos);
        uint64_t m = bits_[i] & ~((~uint64_t(0)) << (int)bit_index(pos));
        if (m) {
            return i*64 + bsr(m);
        } else {
            if (i == 0) return npos;
            i -= 1;
            for (;;) {
                if (bits_[i] != 0) return i*64 + bsr(bits_[i]);
                if (i == 0) return npos;
                i--;
            }
        }
    }

    static const size_t npos = static_cast<size_t>(-1);

protected:
    std::vector<uint64_t> bits_;
    size_t _size = 0;
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
