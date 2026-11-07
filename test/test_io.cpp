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

/**
 * Tests for open_input: decompressing .bz2 input (in parallel for large files, which have
 * many blocks), with several streams, trailing data, and corrupt or truncated files.
 */

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

#include "oink/io.hpp"

#ifdef OINK_HAVE_BZIP2
#include <bzlib.h>
#include <vector>
#include "../src/bzip2_parallel.hpp" // internal, to check that the parallel path is taken
#endif

using namespace pg;

static int failures = 0;

static void
check(bool ok, const std::string& name)
{
    std::cout << (ok ? "ok:   " : "FAIL: ") << name << std::endl;
    if (!ok) failures++;
}

#ifdef OINK_HAVE_BZIP2

/** Pseudo-random bytes (hardly compressible, so the compressed file has many blocks). */
static std::string
random_bytes(size_t n, uint64_t seed)
{
    std::string s(n, '\0');
    for (size_t i = 0; i < n; i++) {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        s[i] = (char)(seed >> 56);
    }
    return s;
}

static std::string
compress(const std::string& data, int level)
{
    unsigned int size = (unsigned int)(data.size() + data.size() / 100 + 600);
    std::string out(size, '\0');
    const int ret = BZ2_bzBuffToBuffCompress(&out[0], &size, const_cast<char*>(data.data()), (unsigned int)data.size(), level, 0, 0);
    if (ret != BZ_OK) throw std::runtime_error("cannot compress test data");
    out.resize(size);
    return out;
}

static std::string
write_temp(const std::string& name, const std::string& contents)
{
    const std::string path = std::string("oink_test_io_") + name + ".bz2";
    std::ofstream f(path, std::ios::binary);
    f << contents;
    return path;
}

static std::string
read_back(const std::string& path)
{
    auto in = open_input(path);
    return std::string((std::istreambuf_iterator<char>(*in)), std::istreambuf_iterator<char>());
}

/** Decompress with the parallel decompressor directly; false if it declines. */
static bool
parallel(const std::string& compressed, std::string& out)
{
    const std::vector<unsigned char> data(compressed.begin(), compressed.end());
    std::vector<std::string> chunks;
    if (!bzip2_decompress_parallel(data, chunks, 4)) return false;
    out.clear();
    for (const auto& c : chunks) out += c;
    return true;
}

static bool
throws(const std::string& path)
{
    try {
        read_back(path);
    } catch (std::runtime_error&) {
        return true;
    }
    return false;
}

static void
test_bzip2()
{
    // 2.5 MB at level 1 (100k blocks): about 25 blocks, over 1 MiB compressed
    const std::string a = random_bytes(2500000, 1);
    const std::string b = random_bytes(1300000, 2) + std::string(200000, 'x');
    const std::string za = compress(a, 1);
    const std::string zb = compress(b, 9);
    check(za.size() >= (1 << 20), "test data is large enough for parallel decompression");

    std::string out;
    check(parallel(za, out) and out == a, "parallel: one stream of many blocks");
    check(parallel(za + zb, out) and out == a + b, "parallel: two streams");
    check(!parallel(za + "not bzip2 data", out), "parallel: declines trailing data (sequential then)");

    std::string p = write_temp("single", za);
    check(read_back(p) == a, "one stream of many blocks");
    std::remove(p.c_str());

    p = write_temp("multi", za + zb);
    check(read_back(p) == a + b, "two concatenated streams (levels 1 and 9)");
    std::remove(p.c_str());

    p = write_temp("trailing", za + "not bzip2 data");
    check(read_back(p) == a, "data after the stream is ignored");
    std::remove(p.c_str());

    std::string corrupt = za;
    corrupt[corrupt.size() / 2] ^= 0x55;
    p = write_temp("corrupt", corrupt);
    check(throws(p), "a corrupt block is an error");
    std::remove(p.c_str());

    p = write_temp("truncated", za.substr(0, za.size() / 2));
    check(throws(p), "a truncated stream is an error");
    std::remove(p.c_str());

    const std::string small = "parity 1;\n0 0 0 1;\n1 1 1 0;\n";
    p = write_temp("small", compress(small, 9));
    check(read_back(p) == small, "a small file (decompressed sequentially)");
    std::remove(p.c_str());
}

#endif

int
main()
{
#ifdef OINK_HAVE_BZIP2
    test_bzip2();
#else
    std::cout << "bzip2 support not built, nothing to test" << std::endl;
#endif
    if (failures) {
        std::cout << failures << " test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all tests passed" << std::endl;
    return 0;
}
