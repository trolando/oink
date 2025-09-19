/*
 * Copyright 2025 Tom van Dijk
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

#include "solvers/gpm.hpp"

#include <queue>
#include <unordered_set>

using namespace pg;

struct BlobHash {
    size_t operator()(const std::vector<std::byte>& v) const noexcept {
        // FNV-1a 64-bit
        size_t h = 1469598103934665603ull;
        for (auto b : v) {
            h ^= size_t(std::to_integer<uint8_t>(b));
            h *= 1099511628211ull;
        }
        return h;
    }
};

using Blob = std::vector<std::byte>;


static Blob snap(Measures* pm, int slot) {
    return pm->store(slot);
}

static void restore(Measures* pm, int slot, const Blob& blob) {
    pm->load(slot, blob);
}

static void test_store_load_equivalence(Measures* pm, int slot) {
    Blob k = snap(pm, slot);
    int other = (slot == 0 ? 1 : 0);

    pm->copy_bot(other);
    restore(pm, other, k);

    if (!pm->eq(slot, other)) {
        std::cerr << "\033[31;1mstore/load mismatch\033[m\n";
        std::cerr << "  original: "; pm->stream(std::cerr, slot);  std::cerr << "\n";
        std::cerr << "  restored: "; pm->stream(std::cerr, other); std::cerr << "\n";
        std::abort();
    }

    if (snap(pm, other) != k) {
        std::cerr << "\033[31;1mblob mismatch after round-trip\033[m\n";
        std::abort();
    }
}



static std::vector<Blob> enumerate_states(Measures* pm, int priorities, bool verbose = false) {
    pm->copy_bot(-4);

    std::queue<Blob> q;
    std::unordered_set<Blob, BlobHash> seen;
    std::vector<Blob> order;

    Blob k0 = snap(pm, -4);
    q.push(k0);
    seen.insert(k0);
    order.push_back(k0);

    while (!q.empty()) {
        Blob cur = q.front(); q.pop();
        restore(pm, -4, cur);

        test_store_load_equivalence(pm, -4);

        for (int a = 0; a < priorities; ++a) {
            pm->copy(-4, -1);
            pm->see(-1, a);

            Blob nxt = snap(pm, -1);
            if (seen.insert(nxt).second) {
                q.push(nxt);
                order.push_back(nxt);
                if (verbose) {
                    std::cout << " + state via a=" << a << ": ";
                    pm->stream(std::cout, -4);
                    std::cout << " -> ";
                    pm->stream(std::cout, -1);
                    std::cout << "\n";
                }
            }
        }
    }
    return order;
}


struct MonoReport {
    size_t pairs_checked = 0;
    size_t succ_checks   = 0;
    size_t violations    = 0;
};

static MonoReport check_monotonicity(Measures* pm, int priorities,
                                     const std::vector<Blob>& states,
                                     size_t max_examples_to_print = 3)
{
    MonoReport rep{};
    size_t printed = 0;

    // We'll use slots:
    //   0 := x, 1 := y, -1 := see_a(x), -4 := see_a(y)
    for (size_t i = 0; i < states.size(); ++i) {
        for (size_t j = 0; j < states.size(); ++j) {
            restore(pm, 0, states[i]);
            restore(pm, 1, states[j]);

            int cmp = pm->compare(0, 1);
            if (cmp <= 0) { // only constrained when x <= y
                ++rep.pairs_checked;

                for (int a = 0; a < priorities; ++a) {
                    // x' = see_a(x)
                    pm->copy(0, -1);
                    pm->see(-1, a);

                    // y' = see_a(y)
                    pm->copy(1, -4);
                    pm->see(-4, a);

                    ++rep.succ_checks;

                    // must have x' <= y'
                    int cmp2 = pm->compare(-1, -4);
                    if (cmp2 > 0) {
                        ++rep.violations;
                        if (printed < max_examples_to_print) {
                            std::cout << "\033[31;1mMonotonicity violation\033[m for a=" << a << "\n";
                            std::cout << "  x:        "; pm->stream(std::cout, 0);  std::cout << "\n";
                            std::cout << "  y:        "; pm->stream(std::cout, 1);  std::cout << "\n";
                            std::cout << "  see_a(x): "; pm->stream(std::cout, -1); std::cout << "\n";
                            std::cout << "  see_a(y): "; pm->stream(std::cout, -4); std::cout << "\n";
                            ++printed;
                        }
                    }
                }
            }
        }
    }

    if (rep.violations == 0) {
        std::cout << "Monotonicity holds for all pairs (x<=y) across all priorities 0.."
                  << (priorities-1) << " (" << rep.pairs_checked << " pairs, "
                  << rep.succ_checks << " successor checks).\n";
    } else {
        std::cout << "Total monotonicity violations: " << rep.violations
                  << " (pairs checked: " << rep.pairs_checked
                  << ", successor checks: " << rep.succ_checks << ").\n";
    }
    return rep;
}


int main(int, char**)
{
    bool print = false;

    // We need some game, just so the measures can be initialized
    int priorities = 7;
    int count = 3;
    Game game(priorities*count);
    {
        game.vec_init();
        int c=0;
        for (int p=0; p<priorities; p++) {
            for (int i=0; i<count; i++) {
                game.init_vertex(c++, p, 0);
            }
        }
        game.sort();
        game.vec_finish();
    }

    std::cout << "Test game has " << priorities << " priorities and " << game.nodecount() << " nodes." << std::endl;

    MeasureKind kinds[] = {
        MeasureKind::Small,
        MeasureKind::Ordered,
        MeasureKind::Succinct
    };

    for (auto kind : kinds) {
        const char* kindName = (kind == MeasureKind::Small)    ? "Small" :
                               (kind == MeasureKind::Ordered)  ? "Ordered" :
                               (kind == MeasureKind::Succinct) ? "Succinct" : "?";

        std::cout << "\n=== MeasureKind: " << kindName << " ===\n";

        for (int player = 0; player < 2; ++player) {
            Measures* pm = new_measure(kind, game, player);
            auto states = enumerate_states(pm, priorities, print);

            std::cout << "[P" << player << "] Reachable unique states from ⊥: "
                      << states.size() << "\n";

            check_monotonicity(pm, priorities, states);
            delete pm;
        }
    }

    /*
    // simple test (just output to manually check)
    if (0) {        
        pm1->copy_bot(-1);
        int arr[] = {1,6,7,5,1,4,5,3,2,1,3,6,3,1,3,3,1,2,1};
        std::cout << "initial: ";
        pm1->stream(std::cout, -1);
        std::cout << std::endl;
        for (int a : arr) {
            pm1->copy(-1, 0);
            pm1->see(-1, a); 
            std::cout << "see a " << a << ": ";
            pm1->stream(std::cout, -1);
            std::cout << " " << pm1->compare(-1, 0);
            std::cout << std::endl;
        }           
    }           
    */

    /*
    // another test to check monotonicity on even and odd measures
    if (0) {
        pm0->copy_bot(0);
        pm0->copy_bot(1);
        if (print) {
            pm0->stream(std::cout, 0);
            std::cout << std::endl;
        }
        while (!pm0->is_top(0)) {
            pm0->copy(0, 1);
            pm0->inc(0);
            if (print) pm0->stream(std::cout, 0);
            if (pm0->compare(0, 1) <= 0) {
                std::cout << " \033[31;1mnot increasing!!\033[m";
            }
            if (print) std::cout << std::endl;
        }
        // return -1;
    }
    */

    return 0;
}
