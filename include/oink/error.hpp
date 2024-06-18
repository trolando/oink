/*
 * Copyright 2019 Tom van Dijk
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

#ifndef PG_ERROR_HPP
#define PG_ERROR_HPP

#include <stdexcept>
#include <sstream>

// #define LOGIC_ERROR { printf("\033[1;7mlogic error %s:%d!\033[m\n", __FILE__, __LINE__); raise(SIGABRT); }
#define THROW_ERROR(msg) { throw pg::Error(msg, __FILE__, __LINE__); }
#define LOGIC_ERROR THROW_ERROR("logic error")

namespace pg {

class Error : public std::exception {
protected:
    const char *msg;
    const char *file;
    const unsigned int line;
    std::string thewhat;
public:
    Error(const char *msg, const char *file, const unsigned int line) : msg(msg), file(file), line(line) {
        std::ostringstream o;
        o << msg << " (at " << file << ":" << line << ")";
        thewhat = o.str();
    }
    ~Error() throw() { }
    const char *what() const noexcept {
        return thewhat.c_str();
    }
};

}

#endif
