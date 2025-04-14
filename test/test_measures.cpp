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

using namespace pg;

int main(int argc, char** argv)
{
    auto measure_kind = MeasureKind::Ordered;

    bool print = false;

    // We need some game, just so the measures can be initialized
    int priorities = 8;
    int count = 2;
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

    Measures *pm0 = new_measure(measure_kind, game, 0);
    Measures *pm1 = new_measure(measure_kind, game, 1);

    // simple test (just output to manually check)
    if (1) {        
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

    // another test to check monotonicity on even and odd measures
    if (1) {
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
        return -1;
    }

    if (0) {
        const int d = game.priority(game.nodecount()-1);
        for (int i=0; i<d; i++) {
            pm0->copy_bot(0);
            pm0->copy_bot(1);
            while (!pm0->is_top(0)) {
                pm0->copy(0, -1);
                pm0->see(-1, i);
                pm0->stream(std::cout, 0);
                std::cout << " see " << i << " => ";
                pm0->stream(std::cout, -1);
                std::cout << std::endl;
                if (pm0->compare(-1, 1) < 0) {
                    std::cout << "ERROR!" << std::endl;
                    return -1;
                }
                pm0->copy(-1, 1);
                pm0->copy(0, -3);
                pm0->inc(0);
                if (pm0->compare(0, -3) <= 0) {
                    std::cout << "ERROR2!" << std::endl;
                    return -1;
                }
            }
        }
        for (int i=0; i<d; i++) {
            pm1->copy_bot(0);
            pm1->copy_bot(1);
            while (!pm1->is_top(0)) {
                pm1->copy(0, -1);
                pm1->see(-1, i);
                pm1->stream(std::cout, 0);
                std::cout << " see " << i << " => ";
                pm1->stream(std::cout, -1);
                std::cout << std::endl;
                if (pm1->compare(-1, 1) < 0) {
                    std::cout << "ERROR!" << std::endl;
                    return -1;
                }
                pm1->copy(-1, 1);
                pm1->copy(0, -3);
                pm1->inc(0); // Inc marks as top too quickly
                if (pm1->compare(0, -3) <= 0) {
                    std::cout << "ERROR2!" << std::endl;
                    return -1;
                }
            }
        }
        return 0;
    }

    return 0;
}
