// SPDX-License-Identifier: Apache-2.0
#include "energy_credit.hpp"
#include "ep_arena.hpp"

#include <iterator>
#include <ostream>
#include <stdexcept>

namespace pg {
using ep::Integer;

struct EnergyCredits::Domain
{
    int player;
    Integer top;
    // For each parity, keep its cumulative magnitude at each present priority.
    // A priority's own magnitude is one plus the opposite total strictly below it.
    std::array<std::vector<std::pair<int, Integer>>, 2> prefix;

    Domain(const Game& game, int player) : player(player)
    {
        Integer totals[2]{};
        for (int v = 0; v < game.nodecount(); ++v) {
            const int p = game.priority(v), color = p & 1;
            ep::next_priority_weight(color, totals);
            if (v + 1 == game.nodecount() || game.priority(v + 1) != p)
                prefix[color].emplace_back(p, totals[color]);
        }
        // A simple prefix against an opponent winning strategy costs at most
        // this player's total positive weight. A larger credit denotes Top.
        top = totals[player] + 1;
    }
};

EnergyCredits::EnergyCredits(const Game& game, int player) :
    domain(std::make_shared<Domain>(game, player)), values(std::size_t(game.nodecount()) + 4) {}

Integer EnergyCredits::progress(const Integer& x, int priority) const
{
    if (x == domain->top) return x;
    const int color = priority & 1;
    const auto& opposite = domain->prefix[1 - color];
    const auto it = std::lower_bound(opposite.begin(), opposite.end(), priority,
        [](const auto& entry, int p) { return entry.first < p; });
    Integer magnitude = 1;
    if (it != opposite.begin()) magnitude += std::prev(it)->second;
    if (color == domain->player) return std::min(Integer(x + magnitude), domain->top);
    return x > magnitude ? Integer(x - magnitude) : Integer(0);
}

bool EnergyCredits::raise(int slot, const Integer& candidate)
{
    if (candidate <= value(slot)) return false;
    values[index(slot)] = candidate;
    return true;
}

void EnergyCredits::copy(int from, int to) { values[index(to)] = value(from); }
void EnergyCredits::copy_bot(int to) { values[index(to)] = 0; }
void EnergyCredits::copy_top(int to) { values[index(to)] = domain->top; }
void EnergyCredits::copy_from(const Measures& other, int from, int to)
{
    const auto& source = dynamic_cast<const EnergyCredits&>(other);
    values[index(to)] = source.is_top(from) ? domain->top : std::min(domain->top, source.value(from));
}
void EnergyCredits::see(int slot, int p) { values[index(slot)] = progress(value(slot), p); }
void EnergyCredits::inc(int slot) { if (!is_top(slot)) ++values[index(slot)]; }
int EnergyCredits::compare(int a, int b) const { return value(a).compare(value(b)); }
bool EnergyCredits::eq(int a, int b) const { return value(a) == value(b); }
bool EnergyCredits::is_top(int slot) const { return value(slot) == domain->top; }
void EnergyCredits::stream(std::ostream& out, int slot) const
{
    if (is_top(slot)) out << "Top";
    else out << value(slot);
}
void EnergyCredits::set(const Measures& other)
{
    const auto& source = dynamic_cast<const EnergyCredits&>(other);
    domain = source.domain;
    values = source.values;
}

std::vector<std::byte> EnergyCredits::store(int slot) const
{
    if (is_top(slot)) return {};
    const auto& x = value(slot);
    std::vector<std::byte> bytes(x == 0 ? 1 : boost::multiprecision::msb(x) / 8 + 1);
    boost::multiprecision::export_bits(x, reinterpret_cast<unsigned char*>(bytes.data()), 8, false);
    return bytes;
}

void EnergyCredits::load(int slot, const std::vector<std::byte>& bytes)
{
    if (bytes.empty()) { copy_top(slot); return; }
    Integer x;
    const auto* first = reinterpret_cast<const unsigned char*>(bytes.data());
    boost::multiprecision::import_bits(x, first, first + bytes.size(), 8, false);
    if (x >= domain->top) throw std::invalid_argument("finite energy credit exceeds the measure bound");
    values[index(slot)] = std::move(x);
}

} // namespace pg
