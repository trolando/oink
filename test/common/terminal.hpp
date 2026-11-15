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
 * Terminal support for the interactive test tools (the solver chooser of
 * test_solvers and bench_solvers, and the live table of bench_solvers).
 *
 * The parts that have no I/O, the key decoder and the text helpers, are kept
 * apart from the Terminal class, so the unit tests can run them without a tty.
 * The Terminal class itself only works on POSIX systems; elsewhere it reports
 * that it is not supported, and the tools fall back to plain output.
 */

#ifndef OINK_TEST_TERMINAL_HPP
#define OINK_TEST_TERMINAL_HPP

#include <csignal>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace ui {

/**
 * A key press or a mouse event read from the terminal.
 */
struct Key
{
    enum Type {
        None, Char, Enter, Tab, Backspace, Delete, Escape, ClearLine,
        Up, Down, Left, Right, Home, End, PageUp, PageDown, Mouse
    };

    Type type = None;
    int ch = 0;           // the character, for Char
    int button = 0;       // the SGR button code of a mouse event (0 left, 64/65 wheel)
    int x = 0;            // the column of a mouse event, counted from 1
    int y = 0;            // the row of a mouse event, counted from 1
    bool release = false; // a mouse button release instead of a press

    Key() = default;
    explicit Key(Type type) : type(type) { }

    static Key character(int c)
    {
        Key key(Char);
        key.ch = c;
        return key;
    }

    bool is(int c) const { return type == Char && ch == c; }
};

/**
 * Turns the bytes read from a terminal into keys: plain characters, the escape
 * sequences of the cursor and editing keys, and SGR mouse reports. It has no
 * I/O, so the tests can feed it any byte sequence, also one split over reads.
 */
class KeyDecoder
{
public:
    void feed(const char* data, size_t size) { buffer.append(data, size); }
    void feed(const std::string& data) { buffer += data; }

    /**
     * Return the next complete key, or a key of type None if the input is empty
     * or ends in an incomplete escape sequence. With <flush> (after a pause in
     * the input), an incomplete escape sequence is reported as Escape.
     */
    Key next(bool flush = false);

    /** Whether undecoded input remains, such as the start of an escape sequence. */
    bool waiting() const { return !buffer.empty(); }

private:
    std::string buffer;
};

/** Replace control characters and non-ASCII bytes by '?'. */
std::string clean(const std::string& text);

/** Clean <text> and pad it to <width> columns, or cut it and end it with '~'. */
std::string fit(const std::string& text, int width);

/** Remove the escape sequences from <text>, to inspect rendered lines in tests. */
std::string strip(const std::string& text);

/**
 * The first visible row of a list of <count> rows shown in a viewport of <rows>
 * rows that started at <first>, such that row <selected> is visible. The
 * viewport only moves when the selection leaves it at an edge.
 */
size_t keep_visible(size_t first, size_t selected, size_t count, size_t rows);

/** Whether file descriptor <fd> is a terminal. */
bool is_tty(int fd);

/** Whether TERM is set to "dumb". */
bool dumb_terminal();

/**
 * Whether output to file descriptor <fd> may use colours and other escape
 * sequences: not with <plain>, only to a terminal, and not with TERM=dumb or
 * when NO_COLOR is set.
 */
bool color_allowed(int fd, bool plain);

/** Escape sequences for the styles of the tools. */
namespace style {
extern const char* const reset;
extern const char* const bold;
extern const char* const dim;
extern const char* const reverse;
extern const char* const cyan;
extern const char* const title;    // bold cyan
extern const char* const good;     // bold green
extern const char* const bad;      // bold red
extern const char* const notice;   // yellow
extern const char* const error;    // red
/** A colour of the 256-colour palette; empty for a negative colour. */
std::string color256(int color);
}

/**
 * The interactive terminal: raw keyboard input from stdin, and a full-screen
 * display on stderr (the alternate screen), so stdout can still be redirected.
 * Only on POSIX systems; elsewhere supported() is false and the constructor
 * throws. The destructor restores the terminal.
 */
class Terminal
{
public:
    /** Whether this platform has the interactive terminal. */
    static bool supported();

    /** Whether stdin and stderr are both terminals and TERM is not dumb. */
    static bool available();

    Terminal();
    ~Terminal();
    Terminal(const Terminal&) = delete;
    Terminal& operator=(const Terminal&) = delete;

    /** The size of the terminal, at least 40 columns and 12 rows. */
    void size(int& width, int& height) const;

    /** Wait at most <timeout_ms> milliseconds for a key; None after the timeout. */
    Key read(int timeout_ms);

    /**
     * Whether a key can be read, waiting at most <timeout_ms> milliseconds for
     * input; nothing is consumed.
     */
    bool pending(int timeout_ms = 0);

    /**
     * Draw the screen: one entry per row, each already fitted to the width.
     * The cursor is shown at <cursor_row>, <cursor_column> (from 1), or hidden
     * when <cursor_row> is 0.
     */
    void draw(const std::vector<std::string>& lines, int cursor_row = 0, int cursor_column = 0);

    /** Report mouse clicks and the wheel as Mouse keys, or stop doing so. */
    void set_mouse(bool on);

    /** Leave the full-screen display and restore the terminal settings. */
    void close();

private:
    struct Saved;
    std::unique_ptr<Saved> saved; // the original terminal settings
    KeyDecoder decoder;
    bool mouse = false;

    bool fill(int timeout_ms);
};

/**
 * While it exists, SIGINT and SIGTERM set <*flag> instead of ending the program,
 * and SIGPIPE is ignored. The previous handlers are restored afterwards.
 */
class InterruptGuard
{
public:
    explicit InterruptGuard(volatile std::sig_atomic_t* flag);
    ~InterruptGuard();
    InterruptGuard(const InterruptGuard&) = delete;
    InterruptGuard& operator=(const InterruptGuard&) = delete;

private:
    volatile std::sig_atomic_t* previous_flag;
    void (*previous_int)(int);
    void (*previous_term)(int);
#ifndef _WIN32
    void (*previous_pipe)(int);
#endif
};

}

#endif
