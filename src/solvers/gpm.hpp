/*
 * Copyright 2022-2026 Tom van Dijk
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
 *
 * The Strahler measures are adapted from the STRPM solver by Vanessa Flügel and
 * Guillermo A. Pérez (https://github.com/vfluegel/oink).
 */

#ifndef GPM_HPP
#define GPM_HPP

#include <optional>
#include <vector>

#include "oink/solver.hpp"

namespace pg {

/**
 * On initialization, all measures must be set to \Bot
 */

class Measures
{
public:
    virtual ~Measures() { }

    // has vertices 0..d plus slots -1, -2

    virtual void copy(int from, int to) = 0; // copy <from> to <to>
    virtual void copy_bot(int to) = 0; // load bottom in <to>
    virtual void copy_top(int to) = 0; // load top in <to>
    virtual void copy_from(const Measures& other, int from, int to) = 0;
    virtual void see(int tgt, int priority) = 0; // let <tgt> see a priority
    virtual void inc(int i) = 0; // smallest increment to <tgt>
    virtual int compare(int i, int j) const = 0; // compare <i> and <j>
    virtual bool eq(int i, int j) const = 0;
    virtual void stream(std::ostream &out, int i) const = 0;
    virtual bool is_top(int tgt) const = 0;
    virtual Measures* clone() const = 0;
    virtual void set(const Measures& other) = 0;
    virtual std::vector<std::byte> store(int tgt) const = 0;
    virtual void load(int tgt, const std::vector<std::byte>& blob) = 0;

    Measures& operator=(const Measures& src) { this->set(src); return *this; }
};

enum class MeasureKind {
    Small,
    Ordered,
    Succinct,
    Strahler
};


class Small : public Measures
{
public:
    Small(const Game &game, const int player);
    Small(const Small &other);
    virtual ~Small();

    void copy(int from, int to) override; // copy <from> to <to>
    void copy_bot(int to) override; // load bottom in <to>
    void copy_top(int to) override; // load top in <to>
    void copy_from(const Measures& other, int from, int to) override;
    void see(int tgt, int priority) override; // let <tgt> see a priority
    int compare(int i, int j) const override; // compare <i> and <j>
    bool eq(int i, int j) const override;
    void stream(std::ostream &out, int i) const override;
    bool is_top(int tgt) const override;
    Measures* clone() const override { return new Small(*this); }
    void set(const Measures& other) override;
    void inc(int i) override; // smallest increment to <tgt>
    std::vector<std::byte> store(int tgt) const override;
    void load(int tgt, const std::vector<std::byte>& blob) override;


private:
    int *data;
    int *counts;
    int n;
    int d;
    int l;
    bitset top;
    int player;
    const Game &game;
};


class Ordered : public Measures
{
public:
    Ordered(const Game &game, const int player);
    Ordered(const Ordered &other);
    virtual ~Ordered();

    void copy(int from, int to) override; // copy <from> to <to>
    void copy_bot(int to) override; // load bottom in <to>
    void copy_top(int to) override; // load top in <to>
    void copy_from(const Measures& other, int from, int to) override;
    void see(int tgt, int priority) override; // let <tgt> see a priority
    int compare(int i, int j) const override; // compare <i> and <j>
    bool eq(int i, int j) const override;
    void stream(std::ostream &out, int i) const override;
    bool is_top(int tgt) const override;
    Measures* clone() const override { return new Ordered(*this); }
    void set(const Measures& other) override;
    void inc(int i) override; // smallest increment to <tgt>
    std::vector<std::byte> store(int tgt) const override;
    void load(int tgt, const std::vector<std::byte>& blob) override;

private:
    int *data;
    int *counts;
    int n;
    int d;
    int l;
    bitset top;
    bitset priorities;
    int player;
    int N; // number of vertices of player priority
    int max_pr; // highest priority of the player
    int max_opp_pr; // highest priority of the opponent

    bool up(int tgt, int priority);
    int val(int tgt) const;
    bool bump(int tgt);
};


class Succinct : public Measures
{
public:
    Succinct(const Game &game, const int player);
    Succinct(const Succinct &other);
    virtual ~Succinct();

    void copy(int from, int to) override; // copy <from> to <to>
    void copy_bot(int to) override; // load bottom in <to>
    void copy_top(int to) override; // load top in <to>
    void copy_from(const Measures& other, int from, int to) override;
    void see(int tgt, int priority) override; // let <tgt> see a priority
    int compare(int i, int j) const override; // compare <i> and <j>
    bool eq(int i, int j) const override;
    void stream(std::ostream &out, int i) const override;
    bool is_top(int tgt) const override;
    Measures* clone() const override { return new Succinct(*this); }
    void set(const Measures& other) override;
    void inc(int i) override; // smallest increment to <tgt>
    std::vector<std::byte> store(int tgt) const override;
    void load(int tgt, const std::vector<std::byte>& blob) override;

private:
    int player; // are we doing even or odd?
    int l; // number of bits: ceil_log2(nodecount())
    int h; // depth of the ordered tree: roughly max_prio/2
           // the even tree (max_prio/2)+1: 1=>1, 2=>2, 3=>2, 4=>3, etc
           // the odd tree (max_prio+1)/2 : 1=>1, 2=>1, 3=>2, 4=>2, etc
    bitset bits; // the bits (l bits per entry)
    int *depth; // depth per bit (assuming it fits in a 32 bit integer)
    long N; // nodecount

    void trunc(int target, int pindex);
    void prog(int target, int pindex); // what is h?


/*
    int *data;
    int *counts;
    int n;
    int d;
    int l;
    bitset top;
    bitset priorities;
    int player;
    int N; // number of vertices of player priority
    int max_pr; // highest priority of the player
    int max_opp_pr; // highest priority of the opponent

    bool up(int tgt, int priority);
    int val(int tgt) const;
    bool bump(int tgt);*/
};





/**
 * Measures that are the leaves of a (k,t)-Strahler-universal tree of height h, from
 * L. Daviaud, M. Jurdzinski, K. S. Thejaswini, "The Strahler Number of a Parity Game",
 * ICALP 2020. The representation (a bitstring with a level per bit, non-decreasing), the
 * successor (prog) and the comparison are from the STRPM solver by Vanessa Flügel with
 * Guillermo A. Pérez (https://github.com/vfluegel/oink), adapted to this interface.
 *
 * A label is a sequence of strings, one per level 0..h-2 (level 0 for the highest
 * priorities); a string is empty (epsilon) or a leading bit followed by non-leading bits,
 * and 0x < epsilon < 1x. A label has exactly k-1 non-empty strings and at most t
 * non-leading bits in total (Definition 21 of the paper). Here t = floor(log2 n) and
 * k = min(t+1, h), so that the tree is (n,h)-universal (Remark 19 of the paper) and
 * the whole game can be solved with it. (STRPM tries pairs up to k = min(t+2, max(h0,h1)),
 * which can be one too small for that, but it solves the game in parts.)
 */
class Strahler final : public Measures
{
public:
    Strahler(const Game &game, const int player);
    Strahler(const Strahler &other) = default;
    virtual ~Strahler() = default;

    void copy(int from, int to) override; // copy <from> to <to>
    void copy_bot(int to) override; // load bottom in <to>
    void copy_top(int to) override; // load top in <to>
    void copy_from(const Measures& other, int from, int to) override;
    void see(int tgt, int priority) override; // let <tgt> see a priority
    int compare(int i, int j) const override; // compare <i> and <j>
    bool eq(int i, int j) const override;
    void stream(std::ostream &out, int i) const override;
    bool is_top(int tgt) const override;
    Measures* clone() const override { return new Strahler(*this); }
    void set(const Measures& other) override;
    void inc(int i) override; // smallest increment to <tgt>
    std::vector<std::byte> store(int tgt) const override;
    void load(int tgt, const std::vector<std::byte>& blob) override;

    int get_k() const { return k; }
    int get_t() const { return t; }
    int get_h() const { return h; }

private:
    int player; // are we doing even or odd?
    int k, t, h; // Strahler number, number of (non-leading) bits, height (levels 0..h-2)
    int cap; // the longest label: k-1 leading bits and t non-leading bits
    int stride; // codes per slot: <cap>, then always END
    long N; // nodecount
    // per slot (slots -4..-1 and the vertices) <stride> codes, one per position: bit 0 at
    // level d is d, bit 1 at level d is 2*END-d, and the positions after the label are END;
    // Top is TOP and then END. Then comparing labels is comparing codes lexicographically
    uint32_t END, TOP;
    std::vector<uint32_t> code;

    // work space for prog and trunc
    std::vector<uint8_t> tb;
    std::vector<int> td;

    void load_tmp(int slot);
    void store_tmp(int slot);
    void prog_tmp(int pindex); // tmp := min { m | m >_p tmp }
    void trunc_tmp(int pindex); // tmp := min { m | m =_p tmp }
    int skip_until_next_level(int i);
};


class GPMSolver : public Solver
{
public:
    GPMSolver(Oink& oink, Game& game, MeasureKind kind = MeasureKind::Ordered);
    virtual ~GPMSolver();

    virtual void run();
    virtual bool parseOptions(std::string&);

private:
    bitset G; // remaining unsolved game
    int *str; // stores currently assigned strategy of each vertex
    MeasureKind measure_kind = MeasureKind::Ordered;
    int lifts = 0;
    int lift_attempts = 0;
    uintqueue Queue;

    bool use_shortcuts = true;

    //bitset R0; // remainder (non-Top, in-Game) for 0
    //bitset R1; // remainder (non-Top, in-Game) for 0
    //bitset B0; // bottom for 0
    //bitset B1; // bottom for 1
    bitset Q0; // queue for 0
    bitset Q1; // queue for 1
    bitset QQ; // next queue
    bitset T0; // top for 0
    bitset T1; // top for 1
    bitset L0; // lifted in last iteration for 0
    bitset L1; // lifted in last iteration for 1

    void attractVertices(int pl, int v, bitset &R, bitset &Z);
    bool lift(Measures &pm, const int player, int v);
    bool update(Measures &pm, const int player);
    void solve(Measures &pm, const int player);
    void shortcuts(const int player, Measures &pm0, Measures &pm1);
};


std::optional<MeasureKind> parse_measure_kind(std::string&);
Measures *new_measure(const MeasureKind kind, const Game &game, const int player);
void stream_measure_kind(std::ostream&, MeasureKind);
void stream_measure_kinds(std::ostream&);

}

#endif
