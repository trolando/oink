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

#include <cassert>
#include <utility>

#include "oink/game_builder.hpp"

namespace pg {

GameBuilder::GameBuilder(int vertex_count)
    : vertex_count_(vertex_count),
      priority_(vertex_count, 0),
      owner_(vertex_count),
      outgoing_(vertex_count)
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
GameBuilder::set_owner(int v, int owner)
{
    assert(v >= 0 and v < vertex_count_);
    owner_[v] = owner ? 1 : 0;
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
    outgoing_[from].push_back(to);
    edge_count_++;
}

void
GameBuilder::truncate(int new_count)
{
    assert(new_count >= 0 and new_count <= vertex_count_);
    for (int v = new_count; v < vertex_count_; v++) edge_count_ -= outgoing_[v].size();
    vertex_count_ = new_count;
    priority_.resize(new_count);
    owner_.resize(new_count);
    if ((int)labels_.size() > new_count) labels_.resize(new_count);
    outgoing_.resize(new_count);
}

Game
GameBuilder::build()
{
    if (vertex_count_ == 0) return Game();

    const size_t ne = edge_count_;

    // The Game accumulator constructor takes labels as nullable pointers; point
    // them at our owned strings (empty string means no label). It copies them,
    // so the borrowed pointers only need to be valid during the call. labels_ is
    // grown lazily, so any index beyond its size has no label.
    std::vector<std::string*> label_ptrs(vertex_count_, nullptr);
    for (int v = 0; v < (int)labels_.size(); v++) {
        if (!labels_[v].empty()) label_ptrs[v] = &labels_[v];
    }

    return Game(static_cast<size_t>(vertex_count_), ne, priority_, owner_, outgoing_, label_ptrs);
}

}
