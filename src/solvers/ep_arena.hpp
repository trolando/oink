// SPDX-License-Identifier: Apache-2.0
#ifndef OINK_EP_ARENA_HPP
#define OINK_EP_ARENA_HPP

#include <algorithm>
#include <array>
#include <cstddef>
#include <iosfwd>
#include <vector>

#include "oink/error.hpp"
#include "oink/game.hpp"
#include "ep_search.hpp"

namespace pg::ep {

// In increasing priority order each vertex outweighs all preceding weights
// of the opposite parity. Equal priorities therefore receive equal weights.
inline ep::Integer next_priority_weight(int color, ep::Integer (&totals)[2])
{
    const ep::Integer magnitude = 1 + totals[1 - color];
    totals[color] += magnitude;
    return color ? -magnitude : magnitude;
}

struct Arena
{
    std::vector<int> vertices;
    std::vector<int> owners;
    std::vector<int> colors;
    std::vector<int> edges;
    std::vector<std::size_t> offsets{0};

    Arena(const Game& game, const bitset& disabled)
    {
        for (int v = 0; v < game.nodecount(); ++v) {
            if (!disabled[v]) vertices.push_back(v);
        }
        std::sort(vertices.begin(), vertices.end(), [&](int a, int b) {
            if (game.priority(a) != game.priority(b)) return game.priority(a) > game.priority(b);
            return a < b;
        });

        std::vector<int> rank(game.nodecount(), -1);
        for (int v = 0; v < size(); ++v) rank[vertices[v]] = v;
        for (int v : vertices) {
            owners.push_back(game.owner(v));
            colors.push_back(game.priority(v) & 1);
            for (int to : game.out_edges(v)) {
                if (rank[to] != -1) edges.push_back(rank[to]);
            }
            if (offsets.back() == edges.size()) THROW_ERROR("EP vertex has no active successor");
            offsets.push_back(edges.size());
        }
    }

    int size() const { return int(vertices.size()); }

    span<const int> successors(int v) const
    {
        return {edges.data() + offsets[v], offsets[v + 1] - offsets[v]};
    }

    std::array<std::vector<int>, 2> branching_vertices() const
    {
        std::array<std::vector<int>, 2> choices;
        for (int v = 0; v < size(); ++v) {
            const auto out = successors(v);
            if (std::any_of(out.begin(), out.end(), [&](int to) { return to != out[0]; }))
                choices[owners[v]].push_back(v);
        }
        return choices;
    }

    template<class Value>
    int best_successor(int v, const std::vector<Value>& values) const
    {
        int best = successors(v)[0];
        for (int to : successors(v)) {
            if (owners[v] == 0 ? values[to] > values[best] : values[to] < values[best]) best = to;
        }
        return best;
    }

    template<class Value>
    void greedy_moves(const std::vector<Value>& values, std::vector<int>& moves) const
    {
        for (int v = 0; v < size(); ++v) moves[v] = best_successor(v, values);
    }

    span<const int> policy_edges(int v, int player, const std::vector<int>& moves) const
    {
        return owners[v] == player ? span<const int>(&moves[v], 1) : successors(v);
    }

    // A rank marker avoids clearing seen between searches. -2 excludes a
    // vertex throughout the search (for example a row with a fixed score).
    template<class Edges>
    bool cycle_at_rank(int top, const Edges& out, std::vector<int>& seen, std::vector<int>& stack) const
    {
        stack.clear();
        stack.push_back(top);
        while (!stack.empty()) {
            const int v = stack.back();
            stack.pop_back();
            for (int to : out(v)) {
                if (to == top) return true;
                if (to < top || seen[to] == top || seen[to] == -2) continue;
                seen[to] = top;
                stack.push_back(to);
            }
        }
        return false;
    }

    std::vector<int> supported_dominions(const ep::Prices& weight, const ep::Prices& low,
                                        const ep::Prices& high, const ep::Integer& unit,
                                        std::vector<int>& moves) const;

    // Check every claimed region, including partial certificates with -1 outside.
    // Then search for a cycle through each wrong-colored vertex, keeping only
    // vertices of that rank or lower. Any losing cycle has such a highest vertex.
    bool certifies(const std::vector<int>& winners, const std::vector<int>& moves) const
    {
        for (int v = 0; v < size(); ++v) {
            if (winners[v] < 0) continue;
            for (int to : policy_edges(v, winners[v], moves)) {
                if (winners[to] != winners[v]) return false;
            }
        }

        std::vector<int> seen(size(), -1), pending;
        for (int top = 0; top < size(); ++top) {
            auto out = [&](int v) { return policy_edges(v, winners[top], moves); };
            if (winners[top] >= 0 && winners[top] != colors[top] && cycle_at_rank(top, out, seen, pending)) return false;
        }
        return true;
    }
};

// The upper endpoint minus an Even winner's minimum energy credit is its
// score. An Odd winner's score is the lower endpoint plus its dual credit.
struct EnergyMeasure
{
    const Arena& arena;
    ep::Prices weight;
    ep::Integer lower = 0, upper = 1;

    explicit EnergyMeasure(const Arena& arena);
    std::vector<int> price_directions(const ep::Prices& trial, const ep::Integer& unit,
                                      std::vector<int>& greedy) const;
    ep::Integer lift(int v, const ep::Integer& successor_score) const;
    ep::Integer lift(int v, const ep::Prices& x) const;
    std::vector<int> solution(const ep::Prices& x, std::vector<int>& moves) const;
};

std::vector<int> lift_measures(const EnergyMeasure& game, std::vector<int>& moves, std::ostream* trace);
std::vector<int> improve_policy(const EnergyMeasure& game, std::vector<int>& moves, std::ostream* trace);
std::vector<int> solve_measures(const EnergyMeasure& game, std::vector<int>& moves, std::ostream* trace);
std::vector<int> sign_solution(const EnergyMeasure& game, std::vector<int>& moves, std::ostream* trace);

}

#endif
