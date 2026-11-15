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
 * A small writer for the JSON objects of the JSON Lines reports. Numbers are
 * written with 17 significant digits, so they read back exactly.
 */

#ifndef OINK_BENCH_JSON_HPP
#define OINK_BENCH_JSON_HPP

#include <iomanip>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace bench {

/** <value> as a JSON string, with quotes and escapes. */
inline std::string
json_string(const std::string& value)
{
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        if (c == '"') out << "\\\"";
        else if (c == '\\') out << "\\\\";
        else if (c == '\n') out << "\\n";
        else if (c == '\r') out << "\\r";
        else if (c == '\t') out << "\\t";
        else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << (unsigned)c << std::dec;
        else out << (char)c;
    }
    out << '"';
    return out.str();
}

/** A number for JSON: 17 significant digits; null if it is not finite. */
inline std::string
json_number(double value)
{
    if (value != value || value - value != 0) return "null";
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

/**
 * Builds one JSON object: Json().add("event", "sample").add("round", 3).str()
 */
class Json
{
public:
    Json& add(const char* key, const std::string& value) { return raw(key, json_string(value)); }
    Json& add(const char* key, const char* value) { return raw(key, json_string(value)); }
    Json& add(const char* key, bool value) { return raw(key, value ? "true" : "false"); }
    Json& add(const char* key, double value) { return raw(key, json_number(value)); }

    template<typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
    Json& add(const char* key, T value)
    {
        return raw(key, std::to_string(value));
    }

    template<typename T>
    Json& add(const char* key, const std::optional<T>& value)
    {
        if (!value) return null(key);
        return add(key, *value);
    }

    Json& null(const char* key) { return raw(key, "null"); }

    /** Add <json>, which is already valid JSON (an object, array or literal). */
    Json& raw(const char* key, const std::string& json)
    {
        if (!text.empty()) text += ',';
        text += json_string(key);
        text += ':';
        text += json;
        return *this;
    }

    std::string str() const { return "{" + text + "}"; }

private:
    std::string text;
};

/** A JSON array of the already formatted <items>. */
inline std::string
json_array(const std::vector<std::string>& items)
{
    std::string result = "[";
    for (size_t i = 0; i < items.size(); i++) {
        if (i) result += ',';
        result += items[i];
    }
    return result + "]";
}

}

#endif
