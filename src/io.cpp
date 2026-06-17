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

#include "oink/io.hpp"

#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <streambuf>
#include <vector>

#ifdef OINK_HAVE_ZLIB
#include <zlib.h>
#endif
#ifdef OINK_HAVE_BZIP2
#include <bzlib.h>
#endif
#ifdef OINK_HAVE_LZMA
#include <cstdint>
#include <lzma.h>
#endif

namespace pg {

namespace {

constexpr std::size_t CHUNK = 1 << 16;

bool has_suffix(const std::string& s, const char* suffix)
{
    const std::string suf(suffix);
    return s.size() >= suf.size() &&
           s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

// An input stream that owns the (decompressed) buffer it reads from.
class buffer_istream : public std::istream {
public:
    explicit buffer_istream(std::string data)
        : std::istream(nullptr), data_(std::move(data)), buf_(data_)
    {
        rdbuf(&buf_);
    }

private:
    struct membuf : std::streambuf {
        explicit membuf(std::string& s)
        {
            char* base = s.empty() ? nullptr : &s[0];
            setg(base, base, base + s.size());
        }
    };
    std::string data_;
    membuf buf_;
};

[[noreturn]] [[maybe_unused]] void unsupported(const char* fmt, const char* lib)
{
    throw std::runtime_error(std::string("compressed input (") + fmt +
        ") is not supported: Oink was built without " + lib);
}

#ifdef OINK_HAVE_ZLIB
std::string inflate_gzip(const std::string& filename)
{
    gzFile f = gzopen(filename.c_str(), "rb");
    if (f == nullptr) throw std::runtime_error("cannot open " + filename);
    std::string out;
    std::vector<char> buf(CHUNK);
    int n;
    while ((n = gzread(f, buf.data(), static_cast<unsigned>(buf.size()))) > 0)
        out.append(buf.data(), n);
    if (n < 0) {
        int errnum = 0;
        const std::string msg = gzerror(f, &errnum);
        gzclose(f);
        throw std::runtime_error("gzip error in " + filename + ": " + msg);
    }
    gzclose(f);
    return out;
}
#endif

#ifdef OINK_HAVE_BZIP2
std::string inflate_bzip2(const std::string& filename)
{
    FILE* fp = std::fopen(filename.c_str(), "rb");
    if (fp == nullptr) throw std::runtime_error("cannot open " + filename);
    int bzerr = BZ_OK;
    BZFILE* bz = BZ2_bzReadOpen(&bzerr, fp, 0, 0, nullptr, 0);
    if (bzerr != BZ_OK) {
        BZ2_bzReadClose(&bzerr, bz);
        std::fclose(fp);
        throw std::runtime_error("cannot read bzip2 stream in " + filename);
    }
    std::string out;
    std::vector<char> buf(CHUNK);
    do {
        const int n = BZ2_bzRead(&bzerr, bz, buf.data(),
                                 static_cast<int>(buf.size()));
        if (n > 0) out.append(buf.data(), n);
    } while (bzerr == BZ_OK);
    const bool ok = (bzerr == BZ_STREAM_END);
    BZ2_bzReadClose(&bzerr, bz);
    std::fclose(fp);
    if (!ok) throw std::runtime_error("bzip2 error in " + filename);
    return out;
}
#endif

#ifdef OINK_HAVE_LZMA
std::string read_file_binary(const std::string& filename)
{
    std::ifstream f(filename, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + filename);
    return std::string((std::istreambuf_iterator<char>(f)),
                       std::istreambuf_iterator<char>());
}

std::string inflate_xz(const std::string& filename)
{
    const std::string compressed = read_file_binary(filename);
    lzma_stream strm = LZMA_STREAM_INIT;
    if (lzma_stream_decoder(&strm, UINT64_MAX, LZMA_CONCATENATED) != LZMA_OK)
        throw std::runtime_error("cannot initialize xz decoder");
    strm.next_in = reinterpret_cast<const uint8_t*>(compressed.data());
    strm.avail_in = compressed.size();
    std::string out;
    std::vector<uint8_t> buf(CHUNK);
    lzma_ret ret;
    do {
        strm.next_out = buf.data();
        strm.avail_out = buf.size();
        ret = lzma_code(&strm, LZMA_FINISH);
        out.append(reinterpret_cast<char*>(buf.data()),
                   buf.size() - strm.avail_out);
    } while (ret == LZMA_OK);
    lzma_end(&strm);
    if (ret != LZMA_STREAM_END)
        throw std::runtime_error("xz error in " + filename);
    return out;
}
#endif

} // namespace

std::unique_ptr<std::istream> open_input(const std::string& filename)
{
    if (has_suffix(filename, ".gz")) {
#ifdef OINK_HAVE_ZLIB
        return std::make_unique<buffer_istream>(inflate_gzip(filename));
#else
        unsupported(".gz", "zlib");
#endif
    }
    if (has_suffix(filename, ".bz2")) {
#ifdef OINK_HAVE_BZIP2
        return std::make_unique<buffer_istream>(inflate_bzip2(filename));
#else
        unsupported(".bz2", "libbz2");
#endif
    }
    if (has_suffix(filename, ".xz")) {
#ifdef OINK_HAVE_LZMA
        return std::make_unique<buffer_istream>(inflate_xz(filename));
#else
        unsupported(".xz", "liblzma");
#endif
    }
    auto f = std::make_unique<std::ifstream>(filename, std::ios::binary);
    if (!*f) throw std::runtime_error("cannot open " + filename);
    return f;
}

std::string supported_input_formats()
{
    std::vector<std::string> fmts;
#ifdef OINK_HAVE_ZLIB
    fmts.emplace_back("gzip (.gz)");
#endif
#ifdef OINK_HAVE_BZIP2
    fmts.emplace_back("bzip2 (.bz2)");
#endif
#ifdef OINK_HAVE_LZMA
    fmts.emplace_back("xz (.xz)");
#endif
    if (fmts.empty()) return "none";
    std::string s = fmts.front();
    for (std::size_t i = 1; i < fmts.size(); ++i) s += ", " + fmts[i];
    return s;
}

} // namespace pg
