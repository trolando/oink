/*
 * Copyright 2024 Tom van Dijk
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

#include <algorithm>
#include <cassert>
#include <utility>

#include "oink/game_builder.hpp"

namespace pg {

GameBuilder::GameBuilder(int vertex_count)
    : vertex_count_(vertex_count),
      priority_(vertex_count, 0),
      owner_(vertex_count)
{
    assert(vertex_count >= 0);
    // labels_ is grown lazily on first set_label, so unlabeled games (the large
    // ones) do not pay for a full vector of empty strings.
}

void
GameBuilder::set_priority(int v, int priority)
{
    assert(v >= 0 and v < vertex_count_);
    priority_[v] = priority;
}

void
GameBuilder::set_owner(int v, Player owner)
{
    assert(v >= 0 and v < vertex_count_);
    owner_[v] = player_index(owner);
}

void
GameBuilder::set_label(int v, std::string label)
{
    assert(v >= 0 and v < vertex_count_);
    if ((int)labels_.size() <= v) labels_.resize(v + 1);
    labels_[v] = std::move(label);
}

void
GameBuilder::add_edge(int from, int to)
{
    assert(from >= 0 and from < vertex_count_);
    assert(to >= 0 and to < vertex_count_);
    edge_from_.push_back(from);
    edge_to_.push_back(to);
}

void
GameBuilder::truncate(int new_count)
{
    assert(new_count >= 0 and new_count <= vertex_count_);
    if (new_count < vertex_count_) {
        // drop the edges of the dropped vertices, keeping the order of the others
        size_t j = 0;
        for (size_t i = 0; i < edge_from_.size(); i++) {
            if (edge_from_[i] < new_count) {
                edge_from_[j] = edge_from_[i];
                edge_to_[j] = edge_to_[i];
                j++;
            }
        }
        edge_from_.resize(j);
        edge_to_.resize(j);
    }
    vertex_count_ = new_count;
    priority_.resize(new_count);
    owner_.resize(new_count);
    if ((int)labels_.size() > new_count) labels_.resize(new_count);
}

Game
GameBuilder::build()
{
    if (vertex_count_ == 0) return Game();

    const int nv = vertex_count_;
    const size_t ne = edge_to_.size();

    Game g(nv, static_cast<int>(ne));
    g.n_vertices = nv;
    g.n_edges = ne;
    std::copy(priority_.begin(), priority_.end(), g._priority.begin());
    g._owner = owner_;

    // Group the edges by source vertex, keeping their order (a counting sort): each vertex
    // gets a block of its out-degree plus a -1 terminator.
    for (const int from : edge_from_) g._outcount[from]++;
    size_t pos = 0;
    for (int v = 0; v < nv; v++) {
        g._firstouts[v] = static_cast<int>(pos);
        pos += g._outcount[v];
        g._outedges[pos++] = -1;
    }
    std::vector<int> next(g._firstouts.begin(), g._firstouts.begin() + nv);
    for (size_t i = 0; i < ne; i++) g._outedges[next[edge_from_[i]]++] = edge_to_[i];
    g.e_size = pos;

    // labels_ is grown lazily, so any index beyond its size has no label
    for (int v = 0; v < (int)labels_.size(); v++) g._label[v] = labels_[v];

    g.is_ordered = std::is_sorted(g._priority.begin(), g._priority.begin() + nv);
    return g;
}

}
