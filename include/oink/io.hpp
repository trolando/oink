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

#ifndef OINK_IO_HPP
#define OINK_IO_HPP

#include <istream>
#include <memory>
#include <string>

namespace pg {

/**
 * Open `filename` for reading, transparently decompressing based on the file
 * name extension:
 *
 *   .gz  -> gzip  (requires zlib)
 *   .bz2 -> bzip2 (requires libbz2)
 *   .xz  -> xz    (requires liblzma)
 *
 * Any other extension is read as a plain, uncompressed file.
 *
 * Compression backends are optional and selected at build time. If `filename`
 * has a recognized compressed extension but the corresponding backend was not
 * compiled in, this throws std::runtime_error.
 *
 * Throws std::runtime_error on I/O or decompression errors.
 */
std::unique_ptr<std::istream> open_input(const std::string& filename);

/**
 * Comma-separated list of compressed input formats compiled into this build
 * (e.g. "gzip (.gz), xz (.xz)"), or "none" if no backend was compiled in.
 */
std::string supported_input_formats();

} // namespace pg

#endif // OINK_IO_HPP
