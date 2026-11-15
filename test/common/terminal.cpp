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

#include "terminal.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

#ifdef _WIN32
#include <io.h>
#else
#include <cerrno>
#include <poll.h>
#include <sys/ioctl.h>
#include <termios.h>
#include <unistd.h>
#endif

namespace ui {

namespace style {
const char* const reset = "\033[0m";
const char* const bold = "\033[1m";
const char* const dim = "\033[2m";
const char* const reverse = "\033[7m";
const char* const cyan = "\033[36m";
const char* const title = "\033[1;36m";
const char* const good = "\033[1;32m";
const char* const bad = "\033[1;31m";
const char* const notice = "\033[33m";
const char* const error = "\033[31m";

std::string
color256(int color)
{
    if (color < 0) return "";
    return "\033[38;5;" + std::to_string(color) + "m";
}
}

static const char ESC = '\033';

// The longest escape sequence we wait for; longer garbage is dropped.
static const size_t MAX_SEQUENCE = 64;

/**
 * Index of the final byte of the control sequence that starts at <from>, or
 * npos if it has not arrived yet. Parameter and intermediate bytes are 0x20 to
 * 0x3f; the final byte is 0x40 to 0x7e.
 */
static size_t
final_byte(const std::string& buffer, size_t from)
{
    for (size_t i = from; i < buffer.size(); i++) {
        const unsigned char c = (unsigned char)buffer[i];
        if (c >= 0x40 && c <= 0x7e) return i;
    }
    return std::string::npos;
}

/** The numbers of the parameters "a;b;c" of a control sequence. */
static std::vector<int>
parameters(const std::string& text)
{
    std::vector<int> result;
    std::istringstream in(text);
    std::string item;
    while (std::getline(in, item, ';')) {
        try {
            result.push_back(item.empty() ? 0 : std::stoi(item));
        } catch (const std::exception&) {
            result.push_back(-1);
        }
    }
    return result;
}

Key
KeyDecoder::next(bool flush)
{
    while (!buffer.empty()) {
        const unsigned char c = (unsigned char)buffer[0];
        if (c != (unsigned char)ESC) {
            buffer.erase(0, 1);
            if (c == '\r' || c == '\n') return Key(Key::Enter);
            if (c == '\t') return Key(Key::Tab);
            if (c == 127 || c == 8) return Key(Key::Backspace);
            if (c == 21) return Key(Key::ClearLine); // Ctrl-U
            return Key::character(c);
        }

        // An escape sequence; wait for its next byte unless the input paused.
        if (buffer.size() == 1) {
            if (!flush) return Key();
            buffer.clear();
            return Key(Key::Escape);
        }

        const char kind = buffer[1];
        if (kind == 'O') {
            // SS3: the cursor keys in application mode, ESC O A
            if (buffer.size() < 3) {
                if (!flush) return Key();
                buffer.erase(0, 1);
                return Key(Key::Escape);
            }
            const char code = buffer[2];
            buffer.erase(0, 3);
            if (code == 'A') return Key(Key::Up);
            if (code == 'B') return Key(Key::Down);
            if (code == 'C') return Key(Key::Right);
            if (code == 'D') return Key(Key::Left);
            if (code == 'H') return Key(Key::Home);
            if (code == 'F') return Key(Key::End);
            continue; // a function key we do not use
        }
        if (kind != '[') {
            // Escape followed by another key (such as Alt-x): report the Escape.
            buffer.erase(0, 1);
            return Key(Key::Escape);
        }

        if (buffer.size() >= 3 && buffer[2] == 'M') {
            // A legacy X10 mouse report: ESC [ M and three bytes. Skipped as a unit.
            if (buffer.size() < 6) {
                if (!flush) return Key();
                buffer.clear();
                continue;
            }
            buffer.erase(0, 6);
            continue;
        }

        const size_t end = final_byte(buffer, 2);
        if (end == std::string::npos) {
            if (buffer.size() > MAX_SEQUENCE) {
                buffer.clear();
                continue;
            }
            if (!flush) return Key();
            buffer.erase(0, 1);
            return Key(Key::Escape);
        }
        const char code = buffer[end];
        const std::string params = buffer.substr(2, end - 2);
        buffer.erase(0, end + 1);

        if (!params.empty() && params[0] == '<') {
            // An SGR mouse report: ESC [ < button ; x ; y M (press) or m (release)
            const auto numbers = parameters(params.substr(1));
            if ((code != 'M' && code != 'm') || numbers.size() != 3) continue;
            if (numbers[0] < 0 || numbers[1] < 1 || numbers[2] < 1) continue;
            Key key(Key::Mouse);
            key.button = numbers[0];
            key.x = numbers[1];
            key.y = numbers[2];
            key.release = code == 'm';
            return key;
        }

        if (code == 'A') return Key(Key::Up);
        if (code == 'B') return Key(Key::Down);
        if (code == 'C') return Key(Key::Right);
        if (code == 'D') return Key(Key::Left);
        if (code == 'H') return Key(Key::Home);
        if (code == 'F') return Key(Key::End);
        if (code == '~') {
            const auto numbers = parameters(params);
            const int number = numbers.empty() ? 0 : numbers[0];
            if (number == 1 || number == 7) return Key(Key::Home);
            if (number == 4 || number == 8) return Key(Key::End);
            if (number == 3) return Key(Key::Delete);
            if (number == 5) return Key(Key::PageUp);
            if (number == 6) return Key(Key::PageDown);
        }
        // Any other control sequence is ignored.
    }
    return Key();
}

std::string
clean(const std::string& text)
{
    std::string result = text;
    for (char& c : result) {
        const unsigned char u = (unsigned char)c;
        if (u < 32 || u >= 127) c = '?';
    }
    return result;
}

std::string
fit(const std::string& text, int width)
{
    if (width <= 0) return "";
    std::string result = clean(text);
    const size_t w = (size_t)width;
    if (result.size() > w) {
        if (w < 2) return result.substr(0, w);
        return result.substr(0, w - 1) + "~";
    }
    return result + std::string(w - result.size(), ' ');
}

std::string
strip(const std::string& text)
{
    std::string result;
    for (size_t i = 0; i < text.size(); i++) {
        if (text[i] != ESC) {
            result += text[i];
            continue;
        }
        if (i + 1 < text.size() && text[i+1] == '[') {
            const size_t end = final_byte(text, i + 2);
            if (end == std::string::npos) break;
            i = end;
        } else {
            i++; // a two-byte escape
        }
    }
    return result;
}

size_t
keep_visible(size_t first, size_t selected, size_t count, size_t rows)
{
    if (count == 0) return 0;
    if (rows < 1) rows = 1;
    if (selected >= count) selected = count - 1;
    // never leave empty rows below the end of the list
    const size_t last_first = count > rows ? count - rows : 0;
    if (first > last_first) first = last_first;
    if (selected < first) return selected;
    if (selected - first >= rows) return selected - rows + 1;
    return first;
}

bool
is_tty(int fd)
{
#ifdef _WIN32
    return _isatty(fd) != 0;
#else
    return isatty(fd) != 0;
#endif
}

bool
dumb_terminal()
{
    const char* term = std::getenv("TERM");
    return term != nullptr && std::string(term) == "dumb";
}

bool
color_allowed(int fd, bool plain)
{
    if (plain || !is_tty(fd) || dumb_terminal()) return false;
    const char* no_color = std::getenv("NO_COLOR");
    return no_color == nullptr || no_color[0] == '\0';
}

/*------------------------------------------------------------------------*/

#ifndef _WIN32

struct Terminal::Saved
{
    termios settings;
};

bool
Terminal::supported()
{
    return true;
}

bool
Terminal::available()
{
    return is_tty(STDIN_FILENO) && is_tty(STDERR_FILENO) && !dumb_terminal();
}

Terminal::Terminal()
{
    if (!available()) throw std::runtime_error("the interactive screen needs a terminal on stdin and stderr");
    auto state = std::make_unique<Saved>();
    if (tcgetattr(STDIN_FILENO, &state->settings) != 0) throw std::runtime_error("cannot read the terminal settings");
    termios raw = state->settings;
    // No line buffering and no echo; Ctrl-C still sends SIGINT.
    raw.c_lflag &= ~(ICANON | ECHO);
    raw.c_cc[VMIN] = 0;
    raw.c_cc[VTIME] = 0;
    if (tcsetattr(STDIN_FILENO, TCSANOW, &raw) != 0) throw std::runtime_error("cannot configure the terminal");
    saved = std::move(state);
    // Alternate screen, hidden cursor, cleared screen.
    std::cerr << "\033[?1049h\033[?25l\033[2J" << std::flush;
}

Terminal::~Terminal()
{
    close();
}

void
Terminal::close()
{
    if (!saved) return;
    set_mouse(false);
    std::cerr << "\033[0m\033[?25h\033[?1049l" << std::flush;
    tcsetattr(STDIN_FILENO, TCSANOW, &saved->settings);
    saved.reset();
}

void
Terminal::size(int& width, int& height) const
{
    width = 100;
    height = 30;
    winsize ws{};
    if (ioctl(STDERR_FILENO, TIOCGWINSZ, &ws) == 0) {
        if (ws.ws_col) width = ws.ws_col;
        if (ws.ws_row) height = ws.ws_row;
    }
    if (width < 40) width = 40;
    if (height < 12) height = 12;
}

/**
 * Wait at most <timeout_ms> for input and add it to the decoder; false if
 * nothing arrived.
 */
bool
Terminal::fill(int timeout_ms)
{
    pollfd fd{STDIN_FILENO, POLLIN, 0};
    const int ready = poll(&fd, 1, timeout_ms);
    if (ready <= 0) return false; // timeout, or interrupted by a signal
    char data[256];
    const ssize_t n = ::read(STDIN_FILENO, data, sizeof(data));
    if (n <= 0) return false;
    decoder.feed(data, (size_t)n);
    return true;
}

Key
Terminal::read(int timeout_ms)
{
    Key key = decoder.next();
    if (key.type != Key::None) return key;
    if (!decoder.waiting() && !fill(timeout_ms)) return Key();
    key = decoder.next();
    // The rest of an escape sequence follows quickly; a lone Escape does not.
    while (key.type == Key::None && decoder.waiting()) {
        if (!fill(30)) return decoder.next(true);
        key = decoder.next();
    }
    return key;
}

bool
Terminal::pending(int timeout_ms)
{
    if (decoder.waiting()) return true;
    return fill(timeout_ms);
}

void
Terminal::draw(const std::vector<std::string>& lines, int cursor_row, int cursor_column)
{
    int width = 0;
    int height = 0;
    size(width, height);
    std::ostringstream out;
    out << "\033[?25l\033[H";
    for (int row = 0; row < height && row < (int)lines.size(); row++) {
        if (row) out << "\r\n";
        out << "\033[2K" << lines[row] << style::reset;
    }
    out << "\033[J";
    if (cursor_row > 0) out << "\033[" << cursor_row << ';' << cursor_column << "H\033[?25h";
    std::cerr << out.str() << std::flush;
}

void
Terminal::set_mouse(bool on)
{
    if (on == mouse) return;
    // Button presses and the wheel (1000) in SGR coordinates (1006).
    std::cerr << (on ? "\033[?1000h\033[?1006h" : "\033[?1000l\033[?1006l") << std::flush;
    mouse = on;
}

#else

struct Terminal::Saved
{
};

bool
Terminal::supported()
{
    return false;
}

bool
Terminal::available()
{
    return false;
}

Terminal::Terminal()
{
    throw std::runtime_error("the interactive screen needs a POSIX terminal; use the command-line options instead");
}

Terminal::~Terminal()
{
}

void
Terminal::close()
{
}

void
Terminal::size(int& width, int& height) const
{
    width = 100;
    height = 30;
}

bool
Terminal::fill(int)
{
    return false;
}

Key
Terminal::read(int)
{
    return Key();
}

bool
Terminal::pending(int)
{
    return false;
}

void
Terminal::draw(const std::vector<std::string>&, int, int)
{
}

void
Terminal::set_mouse(bool)
{
}

#endif

/*------------------------------------------------------------------------*/

static volatile std::sig_atomic_t* interrupt_flag = nullptr;

static void
interrupted(int)
{
    if (interrupt_flag != nullptr) *interrupt_flag = 1;
}

InterruptGuard::InterruptGuard(volatile std::sig_atomic_t* flag)
{
    previous_flag = interrupt_flag;
    interrupt_flag = flag;
    previous_int = std::signal(SIGINT, interrupted);
    previous_term = std::signal(SIGTERM, interrupted);
#ifndef _WIN32
    previous_pipe = std::signal(SIGPIPE, SIG_IGN);
#endif
}

InterruptGuard::~InterruptGuard()
{
    std::signal(SIGINT, previous_int);
    std::signal(SIGTERM, previous_term);
#ifndef _WIN32
    std::signal(SIGPIPE, previous_pipe);
#endif
    interrupt_flag = previous_flag;
}

}
