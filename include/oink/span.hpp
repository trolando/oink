/*
 * Copyright 2024 Tom van Dijk
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

#ifndef SPAN_HPP
#define SPAN_HPP

#include <cstddef>

namespace pg {

/**
 * A minimal non-owning view over a contiguous range (a stand-in for std::span,
 * which is C++20). It holds a pointer and a length and supports range-based for.
 * It is a view: it does not own the data and must not outlive it.
 */
template <typename T>
class span
{
public:
    constexpr span() noexcept : data_(nullptr), size_(0) {}
    constexpr span(T* data, std::size_t size) noexcept : data_(data), size_(size) {}

    [[nodiscard]] constexpr T* data() const noexcept { return data_; }
    [[nodiscard]] constexpr std::size_t size() const noexcept { return size_; }
    [[nodiscard]] constexpr bool empty() const noexcept { return size_ == 0; }

    [[nodiscard]] constexpr T* begin() const noexcept { return data_; }
    [[nodiscard]] constexpr T* end() const noexcept { return data_ + size_; }

    [[nodiscard]] constexpr T& operator[](std::size_t i) const noexcept { return data_[i]; }

private:
    T* data_;
    std::size_t size_;
};

}

#endif
