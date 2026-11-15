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
 * Unit tests for the terminal support of the test tools: decoding keys and
 * mouse reports (also when they arrive split over several reads), fitting text
 * to a width, and the viewport rule of the scrolling lists.
 */

#include <iostream>
#include <string>

#include "common/terminal.hpp"

using namespace ui;

static int failures = 0;

static void
check(const char* name, bool ok)
{
    if (!ok) {
        std::cerr << "FAIL [" << name << "]" << std::endl;
        failures++;
    }
}

/** A short name of a key type, for compact checks. */
static std::string
name(Key::Type type)
{
    switch (type) {
        case Key::Enter: return "<enter>";
        case Key::Tab: return "<tab>";
        case Key::Backspace: return "<bs>";
        case Key::Delete: return "<del>";
        case Key::Escape: return "<esc>";
        case Key::ClearLine: return "<clear>";
        case Key::Up: return "<up>";
        case Key::Down: return "<down>";
        case Key::Left: return "<left>";
        case Key::Right: return "<right>";
        case Key::Home: return "<home>";
        case Key::End: return "<end>";
        case Key::PageUp: return "<pgup>";
        case Key::PageDown: return "<pgdn>";
        case Key::Mouse: return "<mouse>";
        default: return "<?>";
    }
}

/** Decode all of <input> and return the keys as text: characters and key names. */
static std::string
decode(const std::string& input, bool flush = false)
{
    KeyDecoder decoder;
    decoder.feed(input);
    std::string result;
    for (;;) {
        const Key key = decoder.next(flush);
        if (key.type == Key::None) break;
        if (key.type == Key::Char) {
            result += (char)key.ch;
        } else {
            result += name(key.type);
        }
    }
    return result;
}

static void
test_keys()
{
    check("plain characters", decode("ab q") == "ab q");
    check("enter, tab, backspace, ctrl-u", decode("\r\n\t\x7f\x08\x15") == "<enter><enter><tab><bs><bs><clear>");
    check("cursor keys", decode("\033[A\033[B\033[C\033[D") == "<up><down><right><left>");
    check("application cursor keys", decode("\033OA\033OB\033OH\033OF") == "<up><down><home><end>");
    check("editing keys", decode("\033[1~\033[3~\033[4~\033[5~\033[6~\033[H\033[F") == "<home><del><end><pgup><pgdn><home><end>");
    check("modified arrow", decode("\033[1;5C") == "<right>");
    check("unknown sequence is ignored", decode("\033[15~x") == "x");
    check("legacy mouse report is skipped", decode("\033[M abc") == "c");

    // A lone Escape is only reported after a pause (flush).
    check("lone escape waits", decode("\033") == "");
    check("lone escape after a pause", decode("\033", true) == "<esc>");
    check("escape before a key", decode("\033x") == "<esc>x");

    // A sequence split over two reads gives one key.
    {
        KeyDecoder decoder;
        decoder.feed("\033[");
        check("split sequence waits", decoder.next().type == Key::None && decoder.waiting());
        decoder.feed("B");
        check("split sequence completes", decoder.next().type == Key::Down);
        check("split sequence consumed", !decoder.waiting());
    }

    // SGR mouse reports, also split over reads; releases are marked.
    {
        KeyDecoder decoder;
        decoder.feed("\033[<0;1");
        check("split mouse waits", decoder.next().type == Key::None);
        decoder.feed("2;4M\033[<0;12;4m\033[<65;3;7M");
        const Key press = decoder.next();
        check("mouse press", press.type == Key::Mouse && press.button == 0 && press.x == 12 && press.y == 4 && !press.release);
        const Key release = decoder.next();
        check("mouse release", release.type == Key::Mouse && release.release);
        const Key wheel = decoder.next();
        check("mouse wheel down", wheel.type == Key::Mouse && wheel.button == 65 && wheel.x == 3 && wheel.y == 7);
        check("mouse reports consumed", decoder.next().type == Key::None && !decoder.waiting());
    }
    check("malformed mouse report is ignored", decode("\033[<0;1;2;3M.") == ".");
}

static void
test_text()
{
    check("fit pads", fit("ab", 4) == "ab  ");
    check("fit cuts with tilde", fit("abcdef", 4) == "abc~");
    check("fit exact", fit("abcd", 4) == "abcd");
    check("fit tiny width", fit("abc", 1) == "a");
    check("fit zero width", fit("abc", 0).empty());
    check("clean replaces control characters", clean("a\tb\033") == "a?b?");
    check("strip removes styles", strip("\033[1;36mOINK\033[0m \033[38;5;82mx\033[0m") == "OINK x");
    check("color256", style::color256(82) == "\033[38;5;82m" && style::color256(-1).empty());
}

static void
test_viewport()
{
    // the viewport only moves at its edges
    check("selection inside stays", keep_visible(5, 7, 100, 10) == 5);
    check("selection below scrolls", keep_visible(5, 15, 100, 10) == 6);
    check("selection above scrolls", keep_visible(5, 3, 100, 10) == 3);
    check("never past the end", keep_visible(95, 99, 100, 10) == 90);
    check("short list starts at 0", keep_visible(4, 2, 3, 10) == 0);
    check("empty list", keep_visible(4, 0, 0, 10) == 0);
    check("selection clamped", keep_visible(0, 200, 100, 10) == 90);
}

int
test_terminal()
{
    test_keys();
    test_text();
    test_viewport();
    if (failures) std::cerr << failures << " terminal test(s) failed" << std::endl;
    return failures ? 1 : 0;
}
