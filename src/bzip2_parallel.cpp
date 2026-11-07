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

#ifdef OINK_HAVE_BZIP2

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <thread>

#include <bzlib.h>

#include "bzip2_parallel.hpp"

namespace pg {

namespace {

/*
 * The bzip2 format: a stream is "BZh" and a level '1'..'9', then the blocks, then an end
 * marker. A block starts with the 48-bit BLOCK_MAGIC and its 32-bit CRC; the end marker is
 * the 48-bit EOS_MAGIC and the 32-bit stream CRC (combining the block CRCs), padded to a
 * byte. Blocks and the end marker start at any bit offset (most significant bit first).
 */
constexpr uint64_t BLOCK_MAGIC = 0x314159265359ULL;
constexpr uint64_t EOS_MAGIC = 0x177245385090ULL;
constexpr uint64_t MASK48 = 0xFFFFFFFFFFFFULL;

/** Read <n> (at most 64) bits at bit offset <pos>; bits past the end read as 0. */
uint64_t
get_bits(const std::vector<unsigned char>& d, uint64_t pos, int n)
{
    uint64_t v = 0;
    for (int i = 0; i < n; i++, pos++) {
        const uint64_t byte = pos >> 3;
        const int bit = byte < d.size() ? (d[byte] >> (7 - (pos & 7))) & 1 : 0;
        v = (v << 1) | (uint64_t)bit;
    }
    return v;
}

struct Marker {
    uint64_t bit; // bit offset of the marker
    bool eos;     // end of stream (else start of block)
};

/**
 * Find the block and end-of-stream markers starting at bit offsets [8*from, 8*to).
 * A marker at bit offset 8*i+s (s < 8) covers all of byte i+1, so a table of the possible
 * values of that byte skips most positions.
 */
void
scan(const std::vector<unsigned char>& d, size_t from, size_t to, std::vector<Marker>& out)
{
    uint8_t table[256] = {0}; // bit s: a marker at shift s is possible
    for (int s = 0; s < 8; s++) {
        table[(BLOCK_MAGIC >> (32 + s)) & 0xFF] |= (uint8_t)(1 << s);
        table[(EOS_MAGIC >> (32 + s)) & 0xFF] |= (uint8_t)(1 << s);
    }
    const size_t size = d.size();
    if (size < 7) return;
    to = std::min(to, size - 6); // a marker needs 6 or 7 bytes
    for (size_t i = from; i < to; i++) {
        const uint8_t shifts = table[d[i+1]];
        if (shifts == 0) continue;
        uint64_t w = 0; // 64 bits from byte i
        for (int k = 0; k < 8; k++) w = (w << 8) | (i+k < size ? d[i+k] : 0);
        for (int s = 0; s < 8; s++) {
            if (!((shifts >> s) & 1)) continue;
            const uint64_t c = (w >> (16 - s)) & MASK48;
            if (c == BLOCK_MAGIC) out.push_back({8*i + s, false});
            else if (c == EOS_MAGIC) out.push_back({8*i + s, true});
        }
    }
}

struct Block {
    uint64_t begin; // bit offset of the block marker
    uint64_t end;   // bit offset of the next marker
    char level;     // the level of its stream ('1'..'9')
    uint32_t crc;   // the block CRC
};

/** Appends bits (most significant first) to a byte vector. */
struct BitWriter {
    std::vector<unsigned char>& buf;
    uint64_t bit;

    void put(uint64_t v, int n)
    {
        for (int i = n-1; i >= 0; i--, bit++) {
            if ((bit & 7) == 0) buf.push_back(0);
            if ((v >> i) & 1) buf.back() |= (unsigned char)(0x80 >> (bit & 7));
        }
    }
};

/**
 * Decompress <block> as a stream of its own: the stream header, the block, and an end
 * marker whose stream CRC is the block CRC (the CRC of a stream with one block).
 */
bool
decompress_block(const std::vector<unsigned char>& d, const Block& block, std::string& out, std::vector<unsigned char>& buf)
{
    const uint64_t nbits = block.end - block.begin;
    buf.clear();
    buf.reserve(4 + nbits/8 + 12);
    buf.push_back('B');
    buf.push_back('Z');
    buf.push_back('h');
    buf.push_back((unsigned char)block.level);

    // copy the whole bytes of the block, shifted to a byte boundary
    const unsigned s = block.begin & 7;
    size_t k = block.begin >> 3;
    const uint64_t nbytes = nbits / 8;
    for (uint64_t j = 0; j < nbytes; j++, k++) {
        unsigned v = (unsigned)d[k] << s;
        if (s) v |= d[k+1] >> (8 - s);
        buf.push_back((unsigned char)v);
    }
    BitWriter bw{buf, 8*(uint64_t)buf.size()};
    const int rest = (int)(nbits & 7);
    bw.put(get_bits(d, block.begin + 8*nbytes, rest), rest);
    bw.put(EOS_MAGIC, 48);
    bw.put(block.crc, 32);

    bz_stream strm;
    std::memset(&strm, 0, sizeof(strm));
    if (BZ2_bzDecompressInit(&strm, 0, 0) != BZ_OK) return false;
    strm.next_in = reinterpret_cast<char*>(buf.data());
    strm.avail_in = (unsigned)buf.size();
    out.clear();
    out.resize((size_t)(block.level - '0') * 100000 + 4096); // the usual block size, grown if needed
    size_t produced = 0;
    int ret;
    for (;;) {
        if (produced == out.size()) out.resize(out.size() * 2);
        strm.next_out = &out[produced];
        strm.avail_out = (unsigned)(out.size() - produced);
        ret = BZ2_bzDecompress(&strm);
        produced = out.size() - strm.avail_out;
        if (ret != BZ_OK) break;
        if (strm.avail_in == 0 and strm.avail_out != 0) break; // input ended without the end marker
    }
    BZ2_bzDecompressEnd(&strm);
    out.resize(produced);
    return ret == BZ_STREAM_END;
}

}

bool
bzip2_decompress_parallel(const std::vector<unsigned char>& d, std::vector<std::string>& out, unsigned threads)
{
    if (threads == 0) threads = 1;
    const size_t size = d.size();

    // find all markers, splitting the data over the threads
    std::vector<std::vector<Marker>> found(threads);
    {
        std::vector<std::thread> pool;
        const size_t part = (size + threads - 1) / threads;
        for (unsigned t = 0; t < threads; t++) {
            const size_t from = std::min(size, t * part), to = std::min(size, from + part);
            pool.emplace_back([&d, &found, t, from, to] { scan(d, from, to, found[t]); });
        }
        for (auto& th : pool) th.join();
    }
    std::vector<Marker> markers;
    for (auto& f : found) markers.insert(markers.end(), f.begin(), f.end());

    // walk the streams: each block runs up to the next marker, and the markers must line up
    std::vector<Block> blocks;
    size_t m = 0, pos = 0;
    while (pos < size) {
        if (size - pos < 4 or d[pos] != 'B' or d[pos+1] != 'Z' or d[pos+2] != 'h' or d[pos+3] < '1' or d[pos+3] > '9') return false;
        const char level = (char)d[pos+3];
        uint64_t bit = 8*(uint64_t)(pos + 4);
        uint32_t combined = 0;
        for (;;) {
            while (m < markers.size() and markers[m].bit < bit) m++; // (inside a stream header)
            if (m == markers.size() or markers[m].bit != bit) return false;
            if (markers[m].eos) {
                if ((uint32_t)get_bits(d, bit + 48, 32) != combined) return false; // stream CRC
                bit += 80;
                m++;
                break;
            }
            if (m+1 == markers.size()) return false;
            const uint32_t crc = (uint32_t)get_bits(d, bit + 48, 32);
            blocks.push_back({bit, markers[m+1].bit, level, crc});
            combined = ((combined << 1) | (combined >> 31)) ^ crc;
            bit = markers[m+1].bit;
            m++;
        }
        pos = (size_t)((bit + 7) / 8); // the next stream starts at a byte
    }

    // decompress the blocks in parallel (each block checks its own CRC)
    out.assign(blocks.size(), std::string());
    std::atomic<size_t> next{0};
    std::atomic<bool> failed{false};
    {
        std::vector<std::thread> pool;
        const unsigned n = (unsigned)std::min<size_t>(threads, blocks.size());
        for (unsigned t = 0; t < n; t++) {
            pool.emplace_back([&] {
                std::vector<unsigned char> buf;
                for (;;) {
                    const size_t i = next.fetch_add(1);
                    if (i >= blocks.size() or failed.load()) return;
                    if (!decompress_block(d, blocks[i], out[i], buf)) failed.store(true);
                }
            });
        }
        for (auto& th : pool) th.join();
    }
    return !failed.load();
}

}

#endif
