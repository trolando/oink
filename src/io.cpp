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
#include "oink/game.hpp"
#include "oink/solution.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

#ifdef OINK_HAVE_ZLIB
#include <zlib.h>
#endif
#ifdef OINK_HAVE_BZIP2
#include <bzlib.h>
#include <thread>
#include "bzip2_parallel.hpp"
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

// An input stream that owns the (decompressed) data it reads from, as a sequence of
// chunks. A chunk is released once it has been read.
class buffer_istream : public std::istream {
public:
    explicit buffer_istream(std::vector<std::string> chunks)
        : std::istream(nullptr), buf_(std::move(chunks))
    {
        rdbuf(&buf_);
    }

    explicit buffer_istream(std::string data)
        : buffer_istream(std::vector<std::string>(1, std::move(data)))
    {
    }

private:
    struct chunkbuf : std::streambuf {
        explicit chunkbuf(std::vector<std::string> c) : chunks(std::move(c)) { }

        int_type underflow() override
        {
            if (next > 0) std::string().swap(chunks[next-1]); // done with the previous chunk
            while (next < chunks.size()) {
                std::string& s = chunks[next++];
                if (s.empty()) continue;
                setg(&s[0], &s[0], &s[0] + s.size());
                return traits_type::to_int_type(s[0]);
            }
            setg(nullptr, nullptr, nullptr);
            return traits_type::eof();
        }

        std::vector<std::string> chunks;
        size_t next = 0; // the chunk after the current one
    };
    chunkbuf buf_;
};

[[maybe_unused]] std::vector<unsigned char> read_file_bytes(const std::string& filename)
{
    std::ifstream f(filename, std::ios::binary | std::ios::ate);
    if (!f) throw std::runtime_error("cannot open " + filename);
    const std::streamoff size = f.tellg();
    std::vector<unsigned char> data(size > 0 ? (size_t)size : 0);
    f.seekg(0);
    if (!data.empty() and !f.read(reinterpret_cast<char*>(data.data()), size)) throw std::runtime_error("cannot read " + filename);
    return data;
}

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
// Decompress all (concatenated) bzip2 streams in <data>; data after the last stream that
// does not start another stream is ignored, as bzip2 does.
std::string inflate_bzip2_sequential(const std::vector<unsigned char>& data, const std::string& filename)
{
    std::string out;
    std::vector<char> buf(CHUNK);
    size_t pos = 0;
    bool any = false;
    while (pos < data.size()) {
        if (any and (data.size() - pos < 3 or data[pos] != 'B' or data[pos+1] != 'Z' or data[pos+2] != 'h')) break;
        bz_stream strm;
        std::memset(&strm, 0, sizeof(strm));
        if (BZ2_bzDecompressInit(&strm, 0, 0) != BZ_OK) throw std::runtime_error("cannot initialize bzip2 decoder");
        size_t fed = pos; // feed the input in pieces (avail_in is 32 bits)
        auto refill = [&]() {
            const size_t n = std::min<size_t>(data.size() - fed, size_t(1) << 30);
            strm.next_in = const_cast<char*>(reinterpret_cast<const char*>(data.data() + fed));
            strm.avail_in = static_cast<unsigned>(n);
            fed += n;
        };
        refill();
        int ret;
        for (;;) {
            strm.next_out = buf.data();
            strm.avail_out = static_cast<unsigned>(buf.size());
            ret = BZ2_bzDecompress(&strm);
            out.append(buf.data(), buf.size() - strm.avail_out);
            if (ret != BZ_OK) break;
            if (strm.avail_in == 0) {
                if (fed < data.size()) refill();
                else if (strm.avail_out != 0) break; // the input ended inside the stream
            }
        }
        pos = fed - strm.avail_in;
        BZ2_bzDecompressEnd(&strm);
        if (ret != BZ_STREAM_END) throw std::runtime_error("bzip2 error in " + filename);
        any = true;
    }
    if (!any) throw std::runtime_error("bzip2 error in " + filename);
    return out;
}

// Decompress a .bz2 file, in parallel (block by block) if it is large enough.
std::vector<std::string> inflate_bzip2(const std::string& filename)
{
    const std::vector<unsigned char> data = read_file_bytes(filename);
    std::vector<std::string> chunks;
    const unsigned threads = std::max(1u, std::thread::hardware_concurrency());
    if (threads > 1 and data.size() >= (size_t(1) << 20) and bzip2_decompress_parallel(data, chunks, threads)) return chunks;
    chunks.assign(1, inflate_bzip2_sequential(data, filename));
    return chunks;
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

void
write_solution(const Game& game, const Solution& solution, std::ostream& out)
{
    out << "paritysol " << solution.solved().count() << ";" << std::endl;
    for (int i=0; i<game.nodecount(); i++) {
        if (!solution.is_solved(i)) continue;
        out << i << " " << (solution.winner(i) ? "1" : "0");
        if (solution.winner(i) == game.owner(i)) {
            // pgsolver allows one successor: the single strategy, or (for a
            // multi-strategy, where it is the -1 sentinel) a representative move
            int str = solution.strategy(i);
            if (str == -1 and solution.has_multi()) str = first_strategy_edge(game, solution, i);
            if (str != -1) out << " " << str;
        }
        out << ";" << std::endl;
    }
}

void
parse_solution(const Game& game, Solution& solution, std::istream& in)
{
    std::string line;
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        std::string token;
        if (!(ss >> token)) continue;          // ignore empty line
        if (token == "paritysol") continue;    // ignore banner

        int ident = std::stoi(token);
        if (ident < 0 || ident >= game.nodecount()) {
            throw std::runtime_error("node index out of bounds");
        }
        if (solution.is_solved(ident)) throw std::runtime_error("node already solved");

        int w;
        if (!(ss >> w)) throw std::runtime_error("missing winner");
        if (w != 0 && w != 1) throw std::runtime_error("invalid winner");

        int str = -1;
        if (w == game.owner(ident)) {
            if (!(ss >> str)) throw std::runtime_error("missing strategy for winning node");
            // validity of the strategy edge is checked by the verifier
        }
        solution.solve(ident, w, str);
    }
}

} // namespace pg
