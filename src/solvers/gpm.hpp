#ifndef GPM_HPP
#define GPM_HPP

#include <optional>

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

    Measures& operator=(const Measures& src) { this->set(src); return *this; }
};

enum class MeasureKind {
    Small,
    Ordered
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


class GPMSolver : public Solver
{
public:
    GPMSolver(Oink& oink, Game& game);
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
