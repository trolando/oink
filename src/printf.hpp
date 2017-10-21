/* 
   Copyright (c) 2015, 2016 Andreas F. Borchert
   All rights reserved.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
   KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
   WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/*
   This header-only C++11 package provides fmt::printf which
   is intended as a type-safe and extensible drop-in replacement
   for std::printf. The principal idea is to replace

      #include <cstdio>
      // ...
      std::printf(...)

   by

      #include "printf.hpp"       or        #include <printf.hpp>
      // ...
      fmt::printf(...)

   and likewise:

      std::fprintf(fp, ...)       by        fmt::printf(out, ...)
      std::snprintf(s, n, ...)    by        fmt::snprintf(s, n, ...)
      std::wprintf(...)           by        fmt::printf(...)

   where out is an ostream, not a FILE*. As fmt::printf is
   based on variadic template constructs of C++11, this
   is possible in a typesafe way. Consequently, it no
   longer matters for fmt::printf whether you use "%f",
   "%lf", or "%Lf" as format. And all operands are supported
   for which an <<-operator exists.

   fmt::printf uses C++ I/O format flags but makes sure
   that the previous state of manipulators and flags of
   the output stream is restored to its original state
   after its invocation. Any previous state is ignored,
   i.e. fmt::printf("%x", val) will print val in hex
   even if std::cout << std::oct has been used before,
   and the previous octal conversion preference will stay
   in effect for <<-operators after the invocation of fmt::printf.

   Please note that the output format of %p is not
   standardized. As the <<-operator for void* may
   differ from the std::printf behaviour for %p, results
   can be different.

   Wide characters, wide strings, and wide streams and
   their mix are supported. However, the type of the
   format string must match that of the stream.

   Restrictions:
      - The combination of hexfloat with a precision (e.g. "%.2a") 
        is not supported by C++11 (see 22.4.2.2.2 in ISO 14882:2011)
	but supported by std::printf (see 7.21.6.1 in ISO 9899:2011).
	As this implementation depends on the C++11 library, it
	appears hard to find a reasonable workaround for this
	diverting behaviour.

   Alternatives:

   This is not the first attempt to provide printf look
   and feel in a type-safe way for C++:

    - In 1994, Cay S. Horstmann published an article about
      extending the iostreams library in C++ Report where he
      proposed a setformat. Example from his paper:

	 cout << "(" << setformat("%8.2f") << x << "," 
	           << setformat("8.2f") << y << ")" << endl;

      Note that he used one setformat per placeholder as
      C++ at that time did not support variadic templates.
      This paper is available at
	 http://horstmann.com/cpp/iostreams.html

    - The Boost Format library has created an approach
      that does not depend on variadic templates. The
      %-operator is used instead:

      std::cout << boost::format("(x, y) = (%4f, %4f)\n" % x % y;

      See http://www.boost.org/doc/libs/1_59_0/libs/format/doc/format.html

    - There is a proposal by Zhihao Yuan for a printf-like interface for the
      C++ streams library:

      std::cout << std::putf("(x, y) = (%4f, %4f)\n", x, y);

      See http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3506.html
      and https://github.com/lichray/formatxx

    - http://codereview.stackexchange.com/questions/63578/printf-like-formatting-for-stdostream-not-exactly-boostformat
*/

#ifndef FMT_PRINTF_HPP
#define FMT_PRINTF_HPP

#if __cplusplus < 201103L
#error This file requires compiler and library support for the \
ISO C++ 2011 standard.
#else

#include <cassert>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cwchar>
#include <cwctype>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <streambuf>
#include <string>
#include <tuple>
#include <type_traits>

namespace fmt {

namespace impl {

/* type trait to recognize char types which can be distinguished
   from regular numerical types, see also
   http://stackoverflow.com/questions/20958262/char16-t-and-char32-t-types-in-c11
*/
template<typename T> struct is_char : public std::false_type {};
template<> struct is_char<char> : public std::true_type {};
template<> struct is_char<wchar_t> : public std::true_type {};
template<> struct is_char<char16_t> : public std::true_type {};
template<> struct is_char<char32_t> : public std::true_type {};

/* printf is expected to return the number of bytes written;
   the following extensions direct all output to the given
   output stream and count all bytes written */
template<typename CharT, typename Traits = std::char_traits<CharT>>
class counting_ostreambuf : public std::basic_streambuf<CharT, Traits> {
   public:
      counting_ostreambuf(std::basic_streambuf<CharT, Traits>& sbuf) :
	 sbuf(sbuf), nbytes(0) {
      }
      std::streamsize get_count() const {
	 return nbytes;
      }
   protected:
      using Base = std::basic_streambuf<CharT, Traits>;
      using char_type = typename Base::char_type;
      using int_type = typename Base::int_type;
      using traits_type = typename Base::traits_type;

      virtual std::streamsize xsputn(const char_type* s,
	    std::streamsize count) {
	 std::streamsize result = sbuf.sputn(s, count);
	 if (result > 0) nbytes += result;
	 return result;
      }
      virtual int_type overflow(int_type ch) {
	 /* modeled after
	    http://stackoverflow.com/questions/10921761/extending-c-ostream */
	 if (ch == traits_type::eof()) {
	    return traits_type::eof();
	 } else {
	    char_type c = traits_type::to_char_type(ch);
	    return xsputn(&c, 1) == 1? ch: traits_type::eof();
	 }
      }
      virtual int sync() {
	 return sbuf.pubsync();
      }
   private:
      std::basic_streambuf<CharT, Traits>& sbuf;
      std::streamsize nbytes;
};

template<typename CharT, typename Traits = std::char_traits<CharT>>
class counting_ostream : public std::basic_ostream<CharT, Traits> {
   public:
      using Base = std::basic_ostream<CharT, Traits>;
      counting_ostream(std::basic_ostream<CharT, Traits>& out) :
	    Base(&sbuf), sbuf(*(out.rdbuf())) {
	 /* inherit locale from our base stream */
	 this->imbue(out.getloc());
      }
      std::streamsize get_count() const {
	 return sbuf.get_count();
      }
   private:
      counting_ostreambuf<CharT, Traits> sbuf;
};

template<typename CharT, typename Traits = std::char_traits<CharT>>
class uppercase_ostreambuf : public std::basic_streambuf<CharT, Traits> {
   public:
      uppercase_ostreambuf(std::basic_streambuf<CharT, Traits>& sbuf) :
	 sbuf(sbuf) {
      }
   protected:
      using Base = std::basic_streambuf<CharT, Traits>;
      using char_type = typename Base::char_type;
      using int_type = typename Base::int_type;
      using traits_type = typename Base::traits_type;

      virtual std::streamsize xsputn(const char_type* s,
	    std::streamsize count) {
	 for (std::streamsize i = 0; i < count; ++i) {
	    char_type ch = s[i];
	    if (std::islower(ch, this->getloc())) {
	       ch = std::toupper(ch, this->getloc());
	    }
	    if (sbuf.sputc(ch) == traits_type::eof()) {
	       return i;
	    }
	 }
	 return count;
      }

      virtual int_type overflow(int_type ch) {
	 if (ch == traits_type::eof()) {
	    return traits_type::eof();
	 } else {
	    char_type c = traits_type::to_char_type(ch);
	    return xsputn(&c, 1) == 1? ch: traits_type::eof();
	 }
      }
      virtual int sync() {
	 return sbuf.pubsync();
      }
   private:
      std::basic_streambuf<CharT, Traits>& sbuf;
};

template<typename CharT, typename Traits = std::char_traits<CharT>>
class uppercase_ostream : public std::basic_ostream<CharT, Traits> {
   public:
      using Base = std::basic_ostream<CharT, Traits>;
      uppercase_ostream(std::basic_ostream<CharT, Traits>& out) :
	    Base(&sbuf), sbuf(*(out.rdbuf())) {
	 this->copyfmt(out);
	 /* inherit locale from our base stream */
	 this->imbue(out.getloc());
      }
   private:
      uppercase_ostreambuf<CharT, Traits> sbuf;
};

/* std::numpunct extension for the format flag '\''
   that explicitly asks for thousands' grouping characters */
struct thousands_grouping : std::numpunct<char> {
   std::string do_grouping() const {
      return "\3";
   }
};

/* std::numpunct extension that suppresses the
   use of grouping characters,
   this is necessary to conform to std::printf behaviour */
struct suppress_grouping : std::numpunct<char> {
   std::string do_grouping() const {
      return "\0";
   }
};

/* RAII object that saves the current formatting state of the stream
   and makes sure that the state is restored on destruction */
template<typename CharT, typename Traits>
struct format_saver {
   format_saver(std::basic_ios<CharT, Traits>& s) :
      s(s), format_keeper(nullptr) {
      format_keeper.copyfmt(s);
   }
   ~format_saver() {
      s.copyfmt(format_keeper);
   }
   std::basic_ios<CharT, Traits>& s;
   std::basic_ios<CharT, Traits> format_keeper;
};

/* reset the entire format state to its default */
template<typename CharT, typename Traits>
inline void reset_format(std::basic_ios<CharT, Traits>& s) {
   std::basic_ios<CharT, Traits> dflt(nullptr);
   s.copyfmt(dflt);
}

/* internal signed integer type which is used for indices and byte counts */
using integer = std::make_signed<std::size_t>::type;

using flagset = unsigned short;
constexpr flagset is_pointer = 1<<0;
constexpr flagset is_charval = 1<<1;
constexpr flagset is_integer = 1<<2;
constexpr flagset toupper = 1<<3; // when std::uppercase won't cut it
constexpr flagset space_flag = 1<<4; // add space, if non-negative
constexpr flagset plus_flag = 1<<5;
constexpr flagset dyn_width = 1<<6;
constexpr flagset precision = 1<<7; // precision was given
constexpr flagset dyn_precision = 1<<8;
constexpr flagset zero_fill = 1<<9;
constexpr flagset minus_flag = 1<<10;
constexpr flagset special_flag = 1<<11;
constexpr flagset grouping_flag = 1<<12;

/* this structure represents a segment of a format string
   up to and including at most one placeholder */
template<typename CharT>
struct format_segment {
   constexpr format_segment() :
      valid(false),
      beginp(nullptr), endp(nullptr), nextp(nullptr),
      fmtflags(), flags(0), base(0), nof_args(0),
      width(0), precision(0),
      width_index(-1), precision_index(-1), value_index(-1),
      conversion(0) {
   }
   /* valid is set to false when format parsing failed */
   bool valid;
   /* stretch of the format string which is to be
      printed in verbatim; if beginp == endp
      nothing is to be printed */
   const CharT* beginp; const CharT* endp;
   /* where to continue parsing */
   const CharT* nextp;
   /* preliminary set of flags that are to be set
      on the output stream */
   std::ios_base::fmtflags fmtflags;
   /* internal flags */
   flagset flags;
   /* base in case of numerical conversions */
   integer base;
   /* number of arguments that are to be consumed,
      this is between 0 and 2 */
   unsigned short int nof_args;
   /* width and precision, if given within the format */
   std::streamsize width;
   std::streamsize precision;
   /* indices of arguments, where required */
   integer width_index;
   integer precision_index;
   integer value_index;
   /* conversion character, i.e. d o x u etc. */
   CharT conversion;
};

/* parse integer value from format string;
   return false in case of overflows */
template<typename CharT, typename T>
bool parse_integer(const CharT*& format, T& val) {
   T v{};
   CharT ch = *format;
   constexpr T maxval = std::numeric_limits<T>::max();
   constexpr T maxval10 = maxval / 10;
   while (ch >= '0' && ch <= '9') {
      T digit = ch - '0';
      if (v > maxval10) return false;
      v *= 10;
      if (v > maxval - digit) return false;
      v += digit;
      ch = *++format;
   }
   val = v;
   return true;
}

/* parse up to one format specification and
   invoke the respective manipulators for out
   and/or set the corresponding flags */
template<typename CharT>
inline format_segment<CharT>
parse_format_segment(const CharT* format, integer arg_index) {
   format_segment<CharT> result;
   if (!format) return result;

   /* skip everything until we encounter a placeholder
      or the end of the format string */
   result.beginp = format;
   CharT ch = *format;
   while (ch && ch != '%') {
      ch = *++format;
   }
   result.endp = format;

   /* end of format string reached? */
   if (!ch) {
      result.valid = true;
      return result;
   }

   ch = *++format;
   if (!ch) return result; /* format ends with '%' */

   /* process %% */
   if (ch == '%') {
      result.valid = true;
      ++result.endp; /* include first '%' */
      result.nextp = format+1;
      return result;
   }

   /* check if we have an argument index */
   if (ch >= '1' && ch <= '9') {
      const CharT* begin = format;
      integer index;
      if (parse_integer(format, index) && *format == '$') {
	 /* accept argument index */
	 result.value_index = index - 1;
	 ch = *++format;
      } else {
	 /* reset parsing */
	 format = begin; ch = *format;
      }
   }

   /* process conversion flags */
   while (ch == '\'' || ch == '-' || ch == '0' || ch == '+' ||
	 ch == ' ' || ch == '#') {
      switch (ch) {
	 case '\'':
	    result.flags |= grouping_flag;
	    break;
	 case '-':
	    result.flags |= minus_flag;
	    result.fmtflags |= std::ios_base::left;
	    break;
	 case '0': result.flags |= zero_fill; break;
	 case '+':
	    result.flags |= plus_flag;
	    result.fmtflags |= std::ios_base::showpos;
	    break;
	 case ' ': result.flags |= space_flag; break;
	 case '#':
	    result.flags |= special_flag;
	    result.fmtflags |= (std::ios_base::showbase |
	       std::ios_base::showpoint);
	    break;
      }
      ch = *++format;
   }

   if ((result.flags & minus_flag) && (result.flags & zero_fill)) {
      /* if the 0 and - flags both appear, the 0 flag is ignored */
      result.flags &= ~zero_fill;
   }
   if ((result.flags & plus_flag) && (result.flags & space_flag)) {
      /* if the ' ' and '+' flags both appear,
	 the <space> flag shall be ignored */
      result.flags &= ~space_flag;
   }
   /* minimum field width */
   std::streamsize width = 0;
   if (ch == '*') {
      result.flags |= dyn_width; ch = *++format;
      if (ch >= '1' && ch <= '9') {
	 integer index;
	 if (!parse_integer(format, index) || *format != '$') return result;
	 ch = *++format;
	 result.width_index = index - 1;
      } else {
	 result.width_index = arg_index + result.nof_args;
      }
      result.nof_args++;
   } else {
      if (!parse_integer(format, width)) return result;
      ch = *format;
      result.width = width;
   }
   /* precision */
   if (ch == '.') {
      result.flags |= precision;
      ch = *++format;
      std::streamsize precision = 0;
      if (ch == '*') {
	 result.flags |= dyn_precision; ch = *++format;
	 if (ch >= '1' && ch <= '9') {
	    integer index;
	    if (!parse_integer(format, index) || *format != '$') return result;
	    ch = *++format;
	    result.precision_index = index - 1;
	 } else {
	    result.precision_index = arg_index + result.nof_args;
	 }
	 result.nof_args++;
      } else {
	 if (ch >= '0' && ch <= '9') {
	    if (!parse_integer(format, precision)) return result;
	    ch = *format;
	 }
	 result.precision = precision;
      }
      if (result.flags & zero_fill) {
	 /* if a precision is specified, the 0 flag is ignored */
	 result.flags &= ~zero_fill;
      }
   }
   /* skip size specification */
   while (ch == 'l' || ch == 'L' || ch == 'h' ||
	 ch == 'j' || ch == 'z' || ch == 't') {
      ch = *++format;
   }
   /* conversion operation */
   result.conversion = ch;
   switch (ch) {
      case 'd':
      case 'i':
      case 'u':
	 result.flags |= is_integer;
	 result.base = 10;
	 break;
      case 'o':
	 result.flags |= is_integer;
	 result.base = 8;
	 break;
      case 'x':
	 result.flags |= is_integer;
	 result.base = 16;
	 break;
      case 'X':
	 result.fmtflags |= std::ios_base::uppercase;
	 result.flags |= is_integer;
	 result.base = 16;
	 break;
      case 'f':
	 result.fmtflags |= std::ios_base::fixed;
	 result.base = 10;
	 break;
      case 'F':
	 result.fmtflags |= (std::ios_base::fixed | std::ios_base::uppercase);
	 result.flags |= toupper;
	 result.base = 10;
	 break;
      case 'e':
	 result.fmtflags |= std::ios_base::scientific;
	 result.base = 10;
	 break;
      case 'E':
	 result.fmtflags |=
	    std::ios_base::scientific | std::ios_base::uppercase;
	 result.flags |= toupper;
	 result.base = 10;
	 break;
      case 'g':
	 /* default behaviour */
	 result.base = 10;
	 break;
      case 'G':
	 result.fmtflags |= std::ios_base::uppercase;
	 result.flags |= toupper;
	 result.base = 10;
	 break;
      case 'a':
	 result.fmtflags |= std::ios_base::scientific | std::ios_base::fixed;
	 result.base = 16;
	 break;
      case 'A':
	 result.fmtflags |= std::ios_base::scientific |
	    std::ios_base::fixed | std::ios_base::uppercase;
	 result.flags |= toupper;
	 result.base = 16;
	 break;
      case 'p':
	 result.base = 16;
	 result.flags |= is_pointer;
	 break;
      case 'C':
	 /* POSIX extension, equivalent to 'lc' */
      case 'c':
	 result.flags |= is_charval;
	 break;
      case 'S':
	 /* POSIX extension, equivalent to 'ls' */
      case 's':
	 /* when boolean values are printed with %s, we get
	    more readable results; idea taken from N3506 */
	 result.fmtflags |= std::ios_base::boolalpha;
	 break;
      case 'n':
	 /* nothing to be done here */
	 break;
      default:
	 return result;
   }
   if ((result.flags & grouping_flag) && (result.base != 10)) {
      /* grouping is just supported for %i, %d, %u, %f, %F,
         %g, and %G, i.e. all cases with base == 10 */
      result.flags &= ~grouping_flag;
   }
   result.valid = true;
   if (result.value_index < 0) {
      result.value_index = arg_index + result.nof_args;
   }
   result.nof_args++;
   ch = *++format;
   if (ch) {
      result.nextp = format;
   }
   return result;
}

/* similar to std::integer_sequence of C++14 */
template<integer... Is> struct seq {
   typedef seq<Is..., sizeof...(Is)> next;
};
template<integer N> struct gen_seq {
   typedef typename gen_seq<N-1>::type::next type;
};
template<> struct gen_seq<0> {
   typedef seq<> type;
};

/* idea taken from
   http://stackoverflow.com/questions/21062864/optimal-way-to-access-stdtuple-element-in-runtime-by-index
*/

/* apply f on the n-th element of a tuple for compile-time n */
template<integer N, typename Tuple, typename Function>
inline auto apply(const Tuple& tuple, Function&& f)
      -> decltype(f(std::get<N>(tuple))) {
   return f(std::get<N>(tuple));
}

/* helper to apply f on the n-th element of a tuple for runtime n */
template<typename Tuple, typename Function, integer... Is>
inline auto apply(const Tuple& tuple, integer index, Function&& f, seq<Is...>)
      -> decltype(f(std::get<0>(tuple))) {
   using apply_t = decltype(&apply<0, Tuple, Function>);
   static const apply_t apply_functions[] = {&apply<Is, Tuple, Function>...};
   return apply_functions[index](tuple, std::forward<Function>(f));
}

/* apply f on the n-th element of a tuple for runtime n */
template<typename Tuple, typename Function>
inline auto apply(const Tuple& tuple, integer index, Function&& f)
      -> decltype(f(std::get<0>(tuple))) {
   return apply(tuple, index, std::forward<Function>(f),
      typename gen_seq<std::tuple_size<Tuple>::value>::type());
}

/* function object class to extract an integer value by index
   from a tuple */
struct get_value_f {
   template<typename Value>
   typename std::enable_if<std::is_integral<Value>::value, integer>::type
   operator()(Value value) {
      return static_cast<integer>(value);
   }
   /* return -1 when the value is not of integral type */
   template<typename Value>
   typename std::enable_if<!std::is_integral<Value>::value, integer>::type
   operator()(Value) {
      return -1;
   }
};

/* extract an integer value by index from a tuple,
   -1 is returned in case of failures */
template<typename Tuple>
inline integer get_value(const Tuple& tuple, integer index) {
   if (index >= 0 &&
	 index < static_cast<integer>(std::tuple_size<Tuple>::value)) {
      return apply(tuple, index, get_value_f());
   } else {
      return -1;
   }
}

/* set offset value in case of %n */
struct set_value_f {
   set_value_f(std::streamsize offset) : offset(offset) {
   }
   integer operator()(int* ptr) {
      *ptr = static_cast<int>(offset);
      return 0;
   }
   template<typename Value>
   integer operator()(Value) {
      return -1;
   }
   std::streamsize offset;
};

template<typename Tuple>
inline integer set_value(const Tuple& tuple, integer index,
      std::streamsize offset) {
   if (index >= 0 &&
	 index < static_cast<integer>(std::tuple_size<Tuple>::value)) {
      return apply(tuple, index, set_value_f(offset));
   } else {
      return -1;
   }
}

/* general formatted output route */
template<typename CharT, typename Traits, typename Value>
inline typename std::enable_if<
      !std::is_integral<
	 typename std::remove_reference<Value>::type>::value &&
      !std::is_floating_point<
	 typename std::remove_reference<Value>::type>::value, bool>::type
print_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, Value&& value) {
   out << value;
   return !!out;
}

/* formatted output of floating point values */
template<typename CharT, typename Traits, typename Value>
inline typename std::enable_if<
      std::is_floating_point<
	 typename std::remove_reference<Value>::type>::value,
      bool>::type
print_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, Value&& value) {
   if ((fseg.flags & zero_fill) && std::isfinite(value)) {
      out << std::setfill(static_cast<CharT>('0'));
      out.setf(std::ios_base::internal);
   }
   if ((fseg.flags & space_flag) && !std::signbit(value)) {
      if (!out.put(' ')) return false;
      if (fseg.width > 0) {
	 out.width(fseg.width-1);
      }
   }
   if (fseg.flags & toupper) {
      /* the default output operators fail to
	 use uppercase characters in some cases */
      impl::uppercase_ostream<CharT, Traits> fpout(out);
      fpout << value;
   } else {
      out << value;
   }
   return !!out;
}

template<typename Value>
inline integer count_digits(Value value, integer base) {
   if (value == 0) {
      return 1;
   } else {
      integer digits = 0;
      while (value != 0) {
	 value /= base; ++digits;
      }
      return digits;
   }
}

/* formatted output of character values (in case of %c)
   where we got a non-char-type numerical value */
template<typename CharT, typename Traits, typename Value>
inline typename std::enable_if<!is_char<Value>::value, bool>::type
print_char_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>&, Value value) {
   out << static_cast<CharT>(value);
   return !!out;
}

/* formatted output of character values (in case of %c)
   without conversion */
template<typename CharT, typename Traits>
inline bool print_char_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, CharT value) {
   out << value;
   return !!out;
}

/* formatted output of character values (in case of %c)
   that needs to be widened */
template<typename CharT, typename Traits>
inline typename std::enable_if<!std::is_same<char, CharT>::value,
      bool>::type
print_char_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, char value) {
   out << out.widen(value);
   return !!out;
}

/* formatted output of character values (in case of %c)
   where we got a non-char-type numerical value */
template<typename CharT, typename Traits, typename Value>
inline typename std::enable_if<
      is_char<Value>::value &&
      !std::is_same<Value, char>::value &&
      !std::is_same<Value, CharT>::value,
      bool>::type
print_char_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, Value value) {
   auto& f = std::use_facet<std::codecvt<Value, CharT, std::mbstate_t>>(
      out.getloc());
   std::mbstate_t state{};
   std::basic_string<CharT> converted(f.max_length(), 0);
   const Value* from_next;
   CharT* to_next;
   auto result = f.out(state,
      /* from */ &value, &value + 1, from_next,
      /* to */ &converted[0], &converted[converted.size()], to_next);
   if (result == std::codecvt_base::ok) {
      converted.resize(to_next - &converted[0]);
      out << converted;
   } else {
      out.setstate(std::ios_base::failbit);
   }
   return !!out;
}

/* formatted output of integral values
   which possibly need to be converted to characters first
   (in case of %c) */
template<typename CharT, typename Traits, typename Value>
inline typename std::enable_if<
      std::is_integral<typename std::remove_reference<Value>::type>::value,
      bool>::type
print_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, Value&& value) {
   integer padding = 0;
   if (fseg.flags & is_charval) {
      print_char_value(out, fseg, value);
   } else if (fseg.flags & is_integer) {
      if (fseg.flags & zero_fill) {
	 out << std::internal << std::setfill(out.widen('0'));
      } else if (fseg.flags & precision) {
	 integer digits = count_digits(value, fseg.base);
	 integer signwidth = (value < 0) ||
	    (fseg.flags & (plus_flag | space_flag));
	 integer extra = signwidth;
	 if (value != 0 && (fseg.flags & special_flag) && fseg.base == 16) {
	    extra += 2; /* '0x' */
	 }
	 if (fseg.flags & grouping_flag) {
	    extra += digits / 3;
	 }
	 if (fseg.precision > digits) {
	    /* padding with 0s required */
	    if (fseg.width > fseg.precision + extra) {
	       /* manual filling is required */
	       if ((out.flags() & std::ios_base::adjustfield) ==
		     std::ios_base::left) {
		  /* padding has to be postponed */
		  padding = fseg.width - fseg.precision - extra;
	       } else {
		  for (int i = 0; i < fseg.width - fseg.precision - extra;
			++i) {
		     out.put(out.widen(' '));
		  }
	       }
	    }
	    out << std::internal << std::setfill(out.widen('0')) <<
	       std::setw(fseg.precision + extra);
	 }
      }
      if ((fseg.flags & space_flag) && value >= 0) {
	 if (!out.put(' ')) return false;
	 auto width = out.width(0);
	 if (width > 0) {
	    out.width(width-1);
	 }
      }
      /* convert character types to a corresponding integer type */
      using integer = decltype(value + 0);
      if (!(out << static_cast<integer>(value))) return false;
      /* print padding now when it is left adjusted */
      for (int i = 0; i < padding; ++i) {
	 out.put(out.widen(' '));
      }
   } else {
      /* neither %c, %d, %o, %x etc. has been given as expected,
	 we proceed with default behaviour */
      out << value;
   }
   return !!out;
}

/* formatted output of CharT strings;
   precision is honoured */
template<typename CharT, typename Traits>
inline bool print_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, const CharT* value) {
   if (fseg.flags & is_pointer) {
      /* %p given: print pointer value */
      out << static_cast<const void*>(value);
   } else {
      if (fseg.flags & precision) {
	 integer precision = fseg.precision;
	 for (integer i = 0; i < precision; ++i) {
	    if (!value[i]) {
	       precision = i; break;
	    }
	 }
	 integer padding = 0;
	 if (fseg.width > precision) {
	    padding = fseg.width - precision;
	 }
	 bool left = (out.flags() & std::ios_base::adjustfield) ==
		  std::ios_base::left;
	 if (!left) {
	    for (integer i = 0; i < padding; ++i) {
	       out.put(out.widen(' '));
	    }
	 }
	 if (precision > 0) {
	    out.write(value, precision);
	 }
	 if (left) {
	    for (integer i = 0; i < padding; ++i) {
	       out.put(out.widen(' '));
	    }
	 }
      } else {
	 out << value;
      }
   }
   return !!out;
}

/* formatted output of std::nullptr_t strings;
   unfortunately we have no output operator for this type in C++11 */
template<typename CharT, typename Traits>
inline bool print_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, std::nullptr_t value) {
   if (fseg.flags & is_pointer) {
      /* %p given: print pointer value */
      out << static_cast<const void*>(value);
      return !!out;
   } else {
      /* fail this if %p is not given */
      return false;
   }
}

/* formatted output of char strings that need to be widened;
   precision is honoured */
template<typename CharT, typename Traits>
inline typename std::enable_if<!std::is_same<CharT, char>::value, bool>::type
print_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, const char* value) {
   if (fseg.flags & is_pointer) {
      /* %p given: print pointer value */
      out << static_cast<const void*>(value);
   } else {
      integer padding = 0;
      integer len = 0;
      bool left = (out.flags() & std::ios_base::adjustfield) ==
	       std::ios_base::left;
      if (fseg.flags & precision) {
	 len = fseg.precision;
	 for (integer i = 0; i < len; ++i) {
	    if (!value[i]) {
	       len = i; break;
	    }
	 }
      } else {
	 while (value[len]) ++len;
      }
      if (fseg.width > len) {
	 padding = fseg.width - len;
      }
      if (!left) {
	 for (integer i = 0; i < padding; ++i) {
	    out.put(out.widen(' '));
	 }
      }
      for (integer i = 0; i < len; ++i) {
	 out.put(out.widen(value[i]));
      }
      if (left) {
	 for (integer i = 0; i < padding; ++i) {
	    out.put(out.widen(' '));
	 }
      }
   }
   return !!out;
}

/* formatted output of strings that need to be converted;
   precision is honoured */
template<typename CharT, typename Traits, typename Value>
inline typename std::enable_if<
   !std::is_same<CharT, Value>::value && is_char<Value>::value,
   bool>::type
print_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, const Value* value) {
   if (fseg.flags & is_pointer) {
      /* %p given: print pointer value */
      out << static_cast<const void*>(value);
   } else {
      integer len = 0;
      if (fseg.flags & precision) {
	 len = fseg.precision;
	 for (integer i = 0; i < len; ++i) {
	    if (!value[i]) {
	       len = i; break;
	    }
	 }
      } else {
	 while (value[len]) ++len;
      }
      auto& f = std::use_facet<std::codecvt<Value, CharT, std::mbstate_t>>(
	 out.getloc());
      std::mbstate_t state{};
      std::basic_string<CharT> converted(len * f.max_length(), 0);
      const Value* from_next;
      CharT* to_next;
      auto result = f.out(state,
	 /* from */ value, value + len, from_next,
	 /* to */ &converted[0], &converted[converted.size()], to_next);
      if (result == std::codecvt_base::ok) {
	 converted.resize(to_next - &converted[0]);
	 out << converted;
      } else {
	 out.setstate(std::ios_base::failbit);
      }
   }
   return !!out;
}

/* formatted output of non-char pointers that
   have possibly a %p conversion */
template<typename CharT, typename Traits, typename Value>
inline typename std::enable_if<!is_char<Value>::value, bool>::type
print_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, const Value* value) {
   if (fseg.flags & is_pointer) {
      /* print the value of the pointer */
      out << static_cast<const void*>(value);
   } else {
      out << value;
   }
   return !!out;
}

/* formatted output of non-const char pointers
   which are delegated to the const char pointer variants */
template<typename CharT, typename Traits, typename Value>
inline typename std::enable_if<is_char<Value>::value, bool>::type
print_value(std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg, Value* value) {
   return print_value(out, fseg, static_cast<const Value*>(value));
}

template<typename CharT, typename Traits>
struct process_value_f {
   process_value_f(std::basic_ostream<CharT, Traits>& out,
	 const format_segment<CharT>& fseg) :
	 out(out), fseg(fseg) {
   }
   template<typename Value>
   bool operator()(Value&& value) {
      return print_value(out, fseg, std::forward<Value>(value));
   }
   std::basic_ostream<CharT, Traits>& out;
   const format_segment<CharT>& fseg;
};

template<typename Tuple, typename CharT, typename Traits>
inline bool process_value(const Tuple& tuple, integer index,
      std::basic_ostream<CharT, Traits>& out,
      const format_segment<CharT>& fseg) {
   if (index >= 0 &&
	 index < static_cast<integer>(std::tuple_size<Tuple>::value)) {
      return apply(tuple, index, process_value_f<CharT, Traits>(out, fseg));
   } else {
      return false;
   }
}

template<typename Value>
inline std::enable_if<std::is_integral<Value>::value &&
   !std::is_const<Value>::value, bool>
set_value(Value* ptr, std::streamsize value) {
   *ptr = value;
   return true;
}

template<typename Value>
inline bool set_value(Value ptr, std::streamsize value) {
   return false;
}

} // namespace impl

template<typename CharT, typename Traits, typename... Values>
inline int printf(std::basic_ostream<CharT, Traits>& out,
      const CharT* format) {
   impl::counting_ostream<CharT, Traits> cout(out);

   impl::format_segment<CharT> fseg;
   while (format) {
      auto fseg = impl::parse_format_segment(format, 0);
      if (!fseg.valid) return -1;
      if (fseg.nof_args > 0) return -1;
      cout.write(fseg.beginp, fseg.endp - fseg.beginp);
      format = fseg.nextp;
   }
   return cout.get_count();
}

template<typename CharT, typename Traits, typename... Values>
inline int printf(std::basic_ostream<CharT, Traits>& out,
      const CharT* format, Values&&... values) {
   impl::counting_ostream<CharT, Traits> cout(out);
   if (cout.getloc() != std::locale::classic()) {
      cout.imbue(std::locale(cout.getloc(), new impl::suppress_grouping()));
   }
   std::tuple<Values&...> tuple(values...);
   impl::integer nof_args = 0;
   while (format) {
      auto fseg = impl::parse_format_segment(format, nof_args);
      if (!fseg.valid) return -1;
      nof_args += fseg.nof_args;
      if (fseg.endp > fseg.beginp) {
	 cout.write(fseg.beginp, fseg.endp - fseg.beginp);
	 if (!cout) return -1;
      }
      if (fseg.value_index >= 0) {
	 if (fseg.conversion == 'n') {
	    if (impl::set_value(tuple, fseg.value_index,
		  cout.get_count()) < 0) {
	       return -1;
	    }
	 } else {
	    if (fseg.width_index >= 0) {
	       fseg.width = impl::get_value(tuple, fseg.width_index);
	    }
	    if (fseg.precision_index >= 0) {
	       fseg.precision = impl::get_value(tuple, fseg.precision_index);
	    }
	    impl::format_saver<CharT, Traits> fsaver(cout);
	    cout.setf(fseg.fmtflags);
	    cout.setf(fseg.base == 8? std::ios_base::oct  :
		      fseg.base == 10? std::ios_base::dec :
		      fseg.base == 16? std::ios_base::hex :
		      std::ios_base::fmtflags(0), std::ios_base::basefield);
	    if (fseg.width > 0) {
	       cout.width(fseg.width);
	    }
	    if ((fseg.flags & impl::precision) && fseg.precision >= 0) {
	       cout.precision(fseg.precision);
	    }
	    if (fseg.flags & impl::grouping_flag) {
	       cout.imbue(std::locale(cout.getloc(),
		  new impl::thousands_grouping()));
	    }
	    if (!process_value(tuple, fseg.value_index, cout, fseg)) {
	       return -1;
	    }
	 }
      }
      format = fseg.nextp;
   }
   return cout.get_count();
}

template<typename... Values>
inline int printf(const char* format, Values&&... values) {
   return printf(std::cout, format, std::forward<Values>(values)...);
}

template<typename... Values>
inline int printf(const wchar_t* format, Values&&... values) {
   return printf(std::wcout, format, std::forward<Values>(values)...);
}

template<typename... Values>
inline int snprintf(char* s, std::size_t n,
      const char* format, Values&&... values) {
   std::ostringstream os;
   int nbytes = printf(os, format, std::forward<Values>(values)...);
   if (nbytes < 0) return nbytes;
   if (n == 0) return nbytes;
   std::string result(os.str());
   if (nbytes + 1 <= n) {
      std::strcpy(s, result.c_str());
      return nbytes;
   } else {
      std::strcpy(s, result.substr(0, n-1).c_str());
      s[n] = 0;
      return n-1;
   }
}

template<typename... Values>
inline int snprintf(wchar_t* s, std::size_t n,
      const wchar_t* format, Values&&... values) {
   std::wostringstream os;
   int nbytes = printf(os, format, std::forward<Values>(values)...);
   if (nbytes < 0) return nbytes;
   if (n == 0) return nbytes;
   std::wstring result(os.str());
   if (nbytes + 1 <= n) {
      std::wcscpy(s, result.c_str());
      return nbytes;
   } else {
      std::wcscpy(s, result.substr(0, n-1).c_str());
      s[n] = 0;
      return n-1;
   }
}

} // namespace fmt

#endif // of #if __cplusplus < 201103L #else ...
#endif // of #ifndef FMT_PRINTF_HPP
