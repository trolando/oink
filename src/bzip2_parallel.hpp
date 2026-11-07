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

#ifndef OINK_BZIP2_PARALLEL_HPP
#define OINK_BZIP2_PARALLEL_HPP

#include <string>
#include <vector>

namespace pg {

/**
 * Decompress bzip2 data (one or more concatenated streams) using <threads> threads.
 *
 * A bzip2 stream consists of independently compressed blocks, each starting with a 48-bit
 * marker (at any bit offset). This finds the markers, wraps every block in a stream of its
 * own and decompresses the blocks in parallel; the output of block i is <out>[i]. The CRC
 * of every block and of every stream is checked.
 *
 * Returns false if the data does not split into blocks as expected (for example when a
 * marker also occurs inside compressed data) or does not decompress; then <out> is
 * unspecified, and the data should be decompressed sequentially instead, which also gives
 * a proper error for corrupt data.
 */
bool bzip2_decompress_parallel(const std::vector<unsigned char>& data, std::vector<std::string>& out, unsigned threads);

}

#endif
