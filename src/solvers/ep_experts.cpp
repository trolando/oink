// SPDX-License-Identifier: Apache-2.0

#include "ep_experts.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <map>
#include <tuple>
#include <utility>

namespace pg::ep {
namespace {

constexpr unsigned tuple_size = 4;
constexpr unsigned max_guesses = tuple_size + 1;
constexpr unsigned bits_per_exponent = 3;
constexpr unsigned sum_scale = 60;
constexpr unsigned investment_bits = 14;
constexpr unsigned learning_bits = 4;
constexpr int denominator = 1 << (investment_bits + learning_bits);
constexpr int retained = (1 << investment_bits) - 1;
constexpr unsigned outside = 1;

constexpr unsigned factorial[] = {1, 1, 2, 6, 24, 120, 720};
constexpr int prefix_coefficients[6][7] = {
    {0, 60, 0, 0, 0, 0, 0},
    {0, 30, 30, 0, 0, 0, 0},
    {0, 10, 30, 20, 0, 0, 0},
    {0, 0, 15, 30, 15, 0, 0},
    {0, -2, 0, 20, 30, 12, 0},
    {0, 0, -5, 0, 25, 30, 10}
};

using Polynomial = std::map<unsigned, Integer>;

unsigned
binomial(unsigned n, unsigned k)
{
    return factorial[n] / (factorial[k] * factorial[n-k]);
}

unsigned
power_key(unsigned variable, unsigned power)
{
    // Each variable's exponent uses a separate field. The total degree is at most six.
    return power << (bits_per_exponent * (variable-1));
}

unsigned
term_power(unsigned key, unsigned variable)
{
    return (key >> (bits_per_exponent * (variable-1))) & ((1 << bits_per_exponent) - 1);
}

std::array<Integer, 7>
powers(const Integer& base)
{
    std::array<Integer, 7> result;
    result[0] = 1;
    for (unsigned i = 1; i < result.size(); ++i) result[i] = result[i-1] * base;
    return result;
}

struct Bounds
{
    // Variable 0 is fixed at zero. Each entry bounds a difference x[i] - x[j].
    // All implied pairwise bounds are maintained, including after adding a cut.
    std::array<std::array<Integer, max_guesses + 1>, max_guesses + 1> upper;
    unsigned guesses;

    Bounds(unsigned guesses, const Integer& radius) : guesses(guesses)
    {
        for (unsigned i = 1; i <= guesses; ++i) {
            upper[i][0] = radius;
            upper[0][i] = radius;
            for (unsigned j = 1; j <= guesses; ++j) {
                if (j != i) upper[i][j] = 2 * radius;
            }
        }
    }

    bool constrain(unsigned from, unsigned to, const Integer& value)
    {
        if (upper[from][to] <= value) return true;
        if (upper[to][from] < -value) return false;
        for (unsigned i = 0; i <= guesses; ++i) {
            Integer prefix = upper[i][from] + value;
            for (unsigned j = 0; j <= guesses; ++j) {
                Integer through = prefix + upper[to][j];
                if (through < upper[i][j]) upper[i][j] = std::move(through);
            }
        }
        return true;
    }
};

struct Endpoint
{
    unsigned variable;
    Integer offset;
};

void
add_prefix(Polynomial& result, unsigned key, unsigned power, const Integer& coefficient,
           const Endpoint& endpoint, const std::array<Integer, 7>& offset_powers)
{
    // prefix_coefficients[p] is 60 times a polynomial F with F(x) - F(x-1) = x^p.
    // This identity also covers intervals containing negative guesses.
    for (unsigned degree = 1; degree <= power+1; ++degree) {
        int prefix_coefficient = prefix_coefficients[power][degree];
        if (prefix_coefficient == 0) continue;
        Integer factor = coefficient * prefix_coefficient;
        unsigned variable_degree = endpoint.variable == 0 ? 0 : degree;
        for (unsigned k = 0; k <= variable_degree; ++k) {
            const auto& offset_power = offset_powers[degree-k];
            if (offset_power == 0) continue;
            unsigned new_key = k == 0 ? key : key + power_key(endpoint.variable, k);
            result[new_key] += factor * binomial(degree, k) * offset_power;
        }
    }
}

Polynomial
sum_over_range(const Polynomial& poly, unsigned variable, const Endpoint& lower, const Endpoint& upper)
{
    Polynomial result;
    auto lower_powers = powers(lower.offset - 1);
    auto upper_powers = powers(upper.offset);
    for (const auto& [key, coefficient] : poly) {
        if (coefficient == 0) continue;
        unsigned power = term_power(key, variable);
        unsigned remaining_key = key - power_key(variable, power);
        add_prefix(result, remaining_key, power, coefficient, upper, upper_powers);
        add_prefix(result, remaining_key, power, -coefficient, lower, lower_powers);
    }
    return result;
}

Integer
sum_polynomial(const Bounds& region, const Polynomial& poly)
{
    unsigned variable = region.guesses;
    if (variable == 1) {
        Integer sum = 0;
        auto lower = powers(-region.upper[0][1] - 1);
        auto upper = powers(region.upper[1][0]);
        for (const auto& [key, coefficient] : poly) {
            if (coefficient == 0) continue;
            unsigned power = term_power(key, variable);
            Integer prefix = 0;
            for (unsigned degree = 1; degree <= power+1; ++degree) {
                int factor = prefix_coefficients[power][degree];
                if (factor != 0) prefix += factor * (upper[degree] - lower[degree]);
            }
            sum += coefficient * prefix;
        }
        return sum / sum_scale;
    }

    // Projecting a closed difference-bound matrix only drops the eliminated
    // variable. The projected bounds already imply that every lower endpoint is
    // at most every upper endpoint. Split on which endpoints are tight, using their
    // variable indices to break ties so the pieces are disjoint.
    Bounds projected = region;
    --projected.guesses;
    std::array<Endpoint, max_guesses> lower, upper;
    for (unsigned i = 0; i < variable; ++i) {
        lower[i] = {i, -region.upper[i][variable]};
        upper[i] = {i, region.upper[variable][i]};
    }

    Integer sum = 0;
    for (unsigned i = 0; i < variable; ++i) {
        Bounds lower_region = projected;
        bool feasible = true;
        for (unsigned k = 0; feasible && k < variable; ++k) {
            if (k != i) feasible = lower_region.constrain(k, i, lower[i].offset - lower[k].offset - (k < i));
        }
        if (!feasible) continue;
        for (unsigned j = 0; j < variable; ++j) {
            Bounds piece = lower_region;
            feasible = true;
            for (unsigned k = 0; feasible && k < variable; ++k) {
                if (k != j) feasible = piece.constrain(j, k, upper[k].offset - upper[j].offset - (k < j));
            }
            if (feasible) sum += sum_polynomial(piece, sum_over_range(poly, variable, lower[i], upper[j]));
        }
    }
    return sum / sum_scale;
}

using Wealth = std::shared_ptr<const Integer>;

struct Pending
{
    Wealth wealth;
    unsigned maximizers;
    std::shared_ptr<const Pending> parent;
};

struct Cell
{
    Bounds bounds;
    Wealth wealth;
    std::shared_ptr<const Pending> pending;
};

Wealth
settle(const Wealth& parent, const Wealth& child, int gain)
{
    if (parent == child && gain == 0) return parent;
    // In units of 1 / initial_wealth, round (6.6) upwards. The common factor
    // 2^(-14 * level) is omitted, so the child's share enters with weight 2^-14.
    Integer total = *parent * (retained * ((1 << learning_bits) + gain));
    total += *child << learning_bits;
    total += denominator - 1;
    total >>= investment_bits + learning_bits;
    return std::make_shared<Integer>(std::move(total));
}

struct Tuple
{
    std::array<std::size_t, tuple_size> slots{};
    std::array<unsigned, tuple_size> positions{};
    unsigned distinct = 0;
    unsigned permutations = factorial[tuple_size];
    Polynomial multiplicity{{0, 1}};
    // Pool 0 has threshold zero. The other pools are in (level, pivot slot) order.
    std::vector<std::vector<Cell>> pools;

    Tuple(const GroupSlots& sorted, const Integer& radius)
    {
        auto next_powers = powers(radius + 1);
        auto previous_powers = powers(radius);
        for (unsigned i = 0; i < tuple_size;) {
            unsigned first = i;
            while (i < tuple_size && sorted[i] == sorted[first]) ++i;
            unsigned repetitions = i - first;
            slots[distinct] = sorted[first];
            positions[distinct] = first;
            ++distinct;
            permutations /= factorial[repetitions];

            // All options at the same signed slot have identical payoff. Keeping
            // only their maximum guess b gives (b+R+1)^k - (b+R)^k ordered choices
            // for k repeated occurrences of that slot.
            Polynomial product;
            for (unsigned degree = 0; degree < repetitions; ++degree) {
                unsigned power = repetitions - degree;
                Integer coefficient = binomial(repetitions, degree) * (next_powers[power] - previous_powers[power]);
                for (const auto& [key, factor] : multiplicity) {
                    product[key + power_key(distinct, degree)] = factor * coefficient;
                }
            }
            multiplicity = std::move(product);
        }
    }
};

struct PriceCut
{
    Integer non_strict;
    Integer strict;
};

struct Options
{
    unsigned count;
    std::array<unsigned, tuple_size + 1> variable;
    std::array<Integer, tuple_size + 1> price;
    std::array<std::array<PriceCut, tuple_size + 1>, tuple_size + 1> cut;

    Options(const Tuple& tuple, const GroupPrices& prices, const Integer& scale,
            bool has_pivot, const Integer& pivot_price) : count(tuple.distinct)
    {
        // Add the pivot guess to every utility, so the outside option and each
        // slot option become affine in a single guess.
        variable[0] = has_pivot ? count + 1 : 0;
        for (unsigned i = 1; i <= count; ++i) {
            variable[i] = i;
            price[i] = prices[tuple.positions[i-1]] - pivot_price;
        }
        for (unsigned i = 0; i <= count; ++i) {
            for (unsigned j = i+1; j <= count; ++j) {
                Integer delta = price[i] - price[j];
                Integer quotient, remainder;
                boost::multiprecision::divide_qr(delta, scale, quotient, remainder);
                if (remainder < 0) --quotient;
                bool equal = remainder == 0;
                cut[i][j] = {quotient, quotient - equal};
                cut[j][i] = {-quotient - !equal, -quotient - 1};
            }
        }
    }

    bool constrain(Bounds& bounds, unsigned lhs, unsigned rhs, bool strict) const
    {
        const auto& limit = cut[lhs][rhs];
        return bounds.constrain(variable[lhs], variable[rhs], strict ? limit.strict : limit.non_strict);
    }

    bool winner_region(Bounds& bounds, unsigned winner) const
    {
        for (unsigned other = 0; other <= count; ++other) {
            if (other == winner) continue;
            if (!constrain(bounds, other, winner, other < winner)) return false;
        }
        return true;
    }

    Polynomial utility(unsigned winner, const Integer& scale, const Polynomial& multiplicity) const
    {
        Polynomial result;
        for (const auto& [key, coefficient] : multiplicity) {
            result[key + power_key(variable[winner], 1)] += scale * coefficient;
            if (variable[0] != 0) result[key + power_key(variable[0], 1)] -= scale * coefficient;
            if (price[winner] != 0) result[key] -= price[winner] * coefficient;
        }
        return result;
    }
};

void
record_maximizers(std::vector<Cell>& refined, const Cell& cell, Bounds bounds, const Options& options,
                  unsigned winner, unsigned next, unsigned maximizers)
{
    if (next > options.count) {
        auto pending = std::make_shared<Pending>(Pending{cell.wealth, maximizers, cell.pending});
        refined.push_back({std::move(bounds), cell.wealth, std::move(pending)});
        return;
    }

    Bounds tie = bounds;
    if (options.constrain(tie, winner, next, false)) {
        record_maximizers(refined, cell, std::move(tie), options, winner, next+1, maximizers | (1u << next));
    }
    if (options.constrain(bounds, next, winner, true)) {
        record_maximizers(refined, cell, std::move(bounds), options, winner, next+1, maximizers);
    }
}

Integer
cell_objective(const Bounds& bounds, const Options& options,
               const std::array<Polynomial, tuple_size + 1>& utilities)
{
    Integer objective = 0;
    for (unsigned winner = 1; winner <= options.count; ++winner) {
        Bounds region = bounds;
        if (options.winner_region(region, winner)) {
            objective += sum_polynomial(region, utilities[winner]);
        }
    }
    return objective;
}

GroupPrices
group_prices(const Prices& prices, const GroupSlots& slots)
{
    GroupPrices selected;
    for (unsigned i = 0; i < tuple_size; ++i) selected[i] = prices[slots[i]];
    return selected;
}

}

struct ExpertPool::Impl
{
    struct PivotClasses {
        std::vector<std::size_t> label;
        std::size_t count = 1;

        explicit PivotClasses(std::size_t slots) : label(slots) {}

        void observe(const Prices& panel)
        {
            if (count == label.size()) return;
            std::map<std::pair<std::size_t, Integer>, std::size_t> next;
            for (std::size_t s = 0; s < panel.size(); ++s) {
                const auto entry = next.emplace(std::make_pair(label[s], panel[s]), next.size()).first;
                label[s] = entry->second;
            }
            count = next.size();
        }
    };
    using PivotTerms = std::vector<std::pair<std::size_t, std::size_t>>;

    const std::size_t slot_count;
    const Integer radius;
    const std::size_t max_depth;
    const Wealth initial;
    std::size_t level = 0;
    std::vector<GroupSlots> groups;
    std::vector<Tuple> tuples;
    std::vector<std::size_t> coordinate_class;
    std::size_t class_count = 1;
    std::vector<PivotClasses> pivot_class;

    Impl(std::size_t coordinates, const Integer& radius, std::size_t max_depth,
         const Integer& initial_wealth) :
        slot_count(2 * coordinates), radius(radius), max_depth(max_depth),
        initial(std::make_shared<Integer>(initial_wealth)), coordinate_class(coordinates)
    {
        assert(coordinates > 0 && radius >= 0 && initial_wealth > 0);
        // A sorted tuple represents all permutations of its slots and guesses.
        // They have identical utilities and gains throughout the recursion.
        for (std::size_t a = 0; a < slot_count; ++a)
            for (std::size_t b = a; b < slot_count; ++b)
                for (std::size_t c = b; c < slot_count; ++c)
                    for (std::size_t d = c; d < slot_count; ++d) {
                        groups.push_back({a, b, c, d});
                        tuples.emplace_back(groups.back(), radius);
                    }
    }

    std::size_t active_pools() const
    {
        return 1 + level * slot_count;
    }

    // Coordinates in one class have had identical prices and feedback at
    // both signed slots in every previous settlement and pending market.
    template<class Value>
    void refine_classes(const std::vector<Value>& slots)
    {
        if (class_count == coordinate_class.size()) return;
        std::map<std::tuple<std::size_t, Value, Value>, std::size_t> labels;
        for (std::size_t v = 0; v < coordinate_class.size(); ++v) {
            const auto key = std::make_tuple(coordinate_class[v], slots[2*v], slots[2*v + 1]);
            const auto label = labels.emplace(key, labels.size()).first;
            coordinate_class[v] = label->second;
        }
        class_count = labels.size();
    }

    Options options(const Tuple& tuple, std::size_t pool, const GroupPrices& prices,
                    const Integer& scale, const History& history) const
    {
        if (pool == 0) return {tuple, prices, scale, false, 0};
        std::size_t pivot = pool - 1;
        return {tuple, prices, scale, true, history[pivot / slot_count][pivot % slot_count]};
    }

    // A pivot label affects an expert only through its ancestor's price at
    // that level. Its independent guessed bid has the same range at all slots.
    // Equal price histories give identical per-expert settlements, even when
    // the ordinary slot feedback differs. Keep different pivot levels separate.
    PivotTerms pivot_terms(const History& history) const
    {
        PivotTerms terms{{0, 1}};
        for (std::size_t l = 0; l < history.size(); ++l) {
            std::map<std::pair<std::size_t, Integer>, std::size_t> seen;
            for (std::size_t s = 0; s < slot_count; ++s) {
                const auto label = l < pivot_class.size() ? pivot_class[l].label[s] : 0;
                const auto entry = seen.emplace(std::make_pair(label, history[l][s]), terms.size());
                if (entry.second) terms.emplace_back(1 + l * slot_count + s, 1);
                else ++terms[entry.first->second].second;
            }
        }
        return terms;
    }

    Integer group_objective(std::size_t group, const GroupPrices& prices,
                            const Integer& scale, const History& history, const PivotTerms* terms = nullptr) const
    {
        assert(group < tuples.size() && history.size() == level && scale > 0);
#ifndef NDEBUG
        for (unsigned i = 1; i < tuple_size; ++i) {
            if (groups[group][i] == groups[group][i-1]) assert(prices[i] == prices[i-1]);
        }
#endif
        const auto& tuple = tuples[group];
        auto pool_value = [&](std::size_t p) -> Integer {
            const auto market = options(tuple, p, prices, scale, history);
            std::array<Polynomial, tuple_size + 1> utilities;
            for (unsigned winner = 1; winner <= market.count; ++winner) {
                utilities[winner] = market.utility(winner, scale, tuple.multiplicity);
            }

            if (p >= tuple.pools.size()) {
                Bounds bounds(tuple.distinct + (p != 0), radius);
                return cell_objective(bounds, market, utilities) * *initial;
            }
            Integer result = 0;
            for (const auto& cell : tuple.pools[p]) {
                result += cell_objective(cell.bounds, market, utilities) * *cell.wealth;
            }
            return result;
        };
        Integer result = 0;
        if (terms) {
            for (const auto& [p, copies] : *terms) result += pool_value(p) * copies;
        } else {
            for (std::size_t p = 0; p < active_pools(); ++p) result += pool_value(p);
        }
        return result * tuple.permutations;
    }

    Prices minimizing_prices(const Integer& bound, const Integer& scale, const History& history) const
    {
        const auto pivots = pivot_terms(history);
        auto cost = [&](std::size_t g, const GroupPrices& p) { return group_objective(g, p, scale, history, &pivots); };
        if (class_count > 1 && class_count == coordinate_class.size()) return market_prices(class_count, groups, bound, cost);

        // Averaging a signed minimizer over the within-class permutations
        // gives a class-constant minimizer by convexity. In a class-constant
        // trial, replacing each tuple slot by its representative preserves its
        // utility, its gain in every tie and its entire rounded wealth history.
        std::vector<std::size_t> first(class_count, slot_count), counts(class_count);
        for (std::size_t v = 0; v < coordinate_class.size(); ++v) {
            const auto c = coordinate_class[v];
            if (counts[c]++ == 0) first[c] = v;
        }
        std::vector<GroupSlots> reduced;
        std::vector<std::size_t> representative;
        std::vector<Integer> multiplicity;
        for (std::size_t g = 0; g < groups.size(); ++g) {
            GroupSlots slots{};
            Integer copies = 1;
            for (std::size_t i = 0; i < slots.size(); ++i) {
                const auto original = groups[g][i];
                const auto c = coordinate_class[original / 2];
                if (first[c] != original / 2) {
                    copies = 0;
                    break;
                }
                slots[i] = 2 * c + original % 2;
                copies *= counts[c];
            }
            if (copies == 0) continue;
            reduced.push_back(slots);
            representative.push_back(g);
            multiplicity.push_back(std::move(copies));
        }
        // Each ordered slot occurrence has a separate choice of original
        // coordinate. Permutations of the four signed classes and all original
        // pivot choices are already counted in the representative objective.
        if (class_count == 1) {
            // As one signed price varies, utilities have slopes -1, 0 or +1.
            // Integral bids therefore give only integral or half-integral
            // breakpoints. Binary search the exact convex objective on that grid.
            History doubled = history;
            for (auto& panel : doubled) {
                for (auto& price : panel) price *= 2;
            }
            auto value = [&](const Integer& q) {
                Integer total = 0;
                for (std::size_t g = 0; g < reduced.size(); ++g) {
                    GroupPrices trial;
                    for (std::size_t i = 0; i < trial.size(); ++i) {
                        trial[i] = reduced[g][i] == 0 ? q : Integer(-q);
                    }
                    total += multiplicity[g] * group_objective(representative[g], trial, 2 * scale, doubled, &pivots);
                }
                return total;
            };
            Integer price = 0;
            if (bound > 0) {
                const Integer origin = value(0);
                int direction = 0;
                if (value(1) < origin) direction = 1;
                else if (value(-1) < origin) direction = -1;
                if (direction) {
                    Integer lo = 1, hi = 2 * bound;
                    while (lo < hi) {
                        const Integer mid = lo + (hi - lo) / 2;
                        if (value(direction * mid) <= value(direction * (mid + 1))) hi = mid;
                        else lo = mid + 1;
                    }
                    price = direction * lo;
                }
            }
            return Prices(coordinate_class.size(), price);
        }
        const Prices prices = market_prices(class_count, reduced, bound,
            [&](std::size_t g, const GroupPrices& p) -> Integer { return cost(representative[g], p) * multiplicity[g]; });
        Prices expanded(coordinate_class.size());
        for (std::size_t v = 0; v < expanded.size(); ++v) expanded[v] = prices[coordinate_class[v]];
        return expanded;
    }

    void open_probe(const Prices& prices, const Integer& scale, const History& history)
    {
        assert(level < max_depth && history.size() == level && prices.size() == slot_count && scale > 0);
        refine_classes(prices);
        while (pivot_class.size() < history.size()) pivot_class.emplace_back(slot_count);
        for (std::size_t l = 0; l < history.size(); ++l) pivot_class[l].observe(history[l]);
        std::size_t active = active_pools();
        for (std::size_t group = 0; group < tuples.size(); ++group) {
            auto& tuple = tuples[group];
            auto selected = group_prices(prices, groups[group]);
            while (tuple.pools.size() < active) {
                bool has_pivot = !tuple.pools.empty();
                // A newly active pivot pool still has its initial weights.
                // Its pending ancestors had only the outside option at that level.
                std::shared_ptr<const Pending> pending;
                for (std::size_t depth = 0; depth < level; ++depth) {
                    pending = std::make_shared<Pending>(Pending{initial, outside, std::move(pending)});
                }
                tuple.pools.push_back({Cell{Bounds(tuple.distinct + has_pivot, radius), initial, std::move(pending)}});
            }

            for (std::size_t p = 0; p < tuple.pools.size(); ++p) {
                auto& cells = tuple.pools[p];
                if (p >= active) {
                    for (auto& cell : cells) {
                        cell.pending = std::make_shared<Pending>(Pending{cell.wealth, outside, std::move(cell.pending)});
                    }
                    continue;
                }
                const auto market = options(tuple, p, selected, scale, history);
                std::vector<Cell> refined;
                for (const auto& cell : cells) {
                    for (unsigned winner = 0; winner <= market.count; ++winner) {
                        Bounds region = cell.bounds;
                        if (market.winner_region(region, winner)) {
                            record_maximizers(refined, cell, std::move(region), market, winner, winner+1, 1u << winner);
                        }
                    }
                }
                cells = std::move(refined);
            }
        }
        ++level;
    }

    void close_probe(const std::vector<int>& feedback)
    {
        assert(level > 0 && feedback.size() == slot_count);
        refine_classes(feedback);
        for (auto& tuple : tuples) {
            for (auto& cells : tuple.pools) {
                for (auto& cell : cells) {
                    const auto& pending = *cell.pending;
                    int gain = (pending.maximizers & outside) ? 0 : 1;
                    for (unsigned i = 0; i < tuple.distinct; ++i) {
                        if (pending.maximizers & (1u << (i+1))) gain = std::min(gain, feedback[tuple.slots[i]]);
                    }
                    cell.wealth = settle(pending.wealth, cell.wealth, gain);
                    cell.pending = pending.parent;
                }
            }
        }
        --level;
    }
};

ExpertPool::ExpertPool(std::size_t coordinates, const Integer& radius, std::size_t max_depth,
                       const Integer& initial_wealth) :
    impl_(std::make_unique<Impl>(coordinates, radius, max_depth, initial_wealth))
{}

ExpertPool::~ExpertPool() = default;

Integer
ExpertPool::count(std::size_t coordinates, const Integer& radius, std::size_t max_depth)
{
    Integer choices = 2 * Integer(coordinates) * (2 * radius + 1);
    Integer pairs = choices * choices;
    return pairs * pairs * (1 + max_depth * choices);
}

const std::vector<GroupSlots>&
ExpertPool::group_slots() const
{
    return impl_->groups;
}

Integer
ExpertPool::group_objective(std::size_t group, const GroupPrices& prices,
                            const Integer& scale, const History& history) const
{
    return impl_->group_objective(group, prices, scale, history);
}

Integer
ExpertPool::objective(const Prices& prices, const Integer& scale, const History& history) const
{
    assert(prices.size() == impl_->slot_count);
    Integer result = 0;
    for (std::size_t group = 0; group < impl_->groups.size(); ++group) {
        result += group_objective(group, group_prices(prices, impl_->groups[group]), scale, history);
    }
    return result;
}

Prices
ExpertPool::minimizing_prices(const Integer& bid_bound, const Integer& scale, const History& history) const
{
    return impl_->minimizing_prices(bid_bound, scale, history);
}

void
ExpertPool::open_probe(const Prices& prices, const Integer& scale, const History& history)
{
    impl_->open_probe(prices, scale, history);
}

void
ExpertPool::close_probe(const std::vector<int>& feedback)
{
    impl_->close_probe(feedback);
}

}
