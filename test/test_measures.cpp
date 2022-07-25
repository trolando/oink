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

namespace pg {

    int main(int argc, char** argv)
    {
        auto measure_kind = MeasureKind::Small;

        Measures *pm0 = new_measure(measure_kind, game, 0);
        Measures *pm1 = new_measure(measure_kind, game, 1);

        // simple test (just output to manually check)
        if (0) {        
            pm1->copy_bot(-1);
            int arr[] = {1,6,7,5,1,4,5,3,2,1,3,2,3,1,3,3,1,2,1};
            for (int a : arr) {
                pm1->see(-1, a); 
                logger << "see a " << a << ": ";
                pm1->stream(logger, -1);
                logger << std::endl;
            }           
        }           
        // another test to check monotonicity on even and odd measures
        if (0) {
            pm0->copy_bot(0);
            pm0->copy_bot(1);
            while (!pm0->is_top(0)) {
                pm0->see(0, 0);
                logger << "see " << 0 << " => ";
                pm0->stream(logger, 0);
                //logger << " val: " << pm0->val(0);
                logger << std::endl;
            }
            return;
        }


        if (0) {
            const int d = priority(nodecount()-1);
            for (int i=0; i<d; i++) {
                pm0->copy_bot(0);
                pm0->copy_bot(1);
                while (!pm0->is_top(0)) {
                    pm0->copy(0, -1);
                    pm0->see(-1, i);
                    pm0->stream(logger, 0);
                    logger << " see " << i << " => ";
                    pm0->stream(logger, -1);
                    logger << std::endl;
                    if (pm0->compare(-1, 1) < 0) {
                        logger << "ERROR!" << std::endl;
                        return;
                    }
                    pm0->copy(-1, 1);
                    pm0->copy(0, -3);
                    pm0->inc(0);
                    if (pm0->compare(0, -3) <= 0) {
                        logger << "ERROR2!" << std::endl;
                        return;
                    }
                }
            }
            for (int i=0; i<d; i++) {
                pm1->copy_bot(0);
                pm1->copy_bot(1);
                while (!pm1->is_top(0)) {
                    pm1->copy(0, -1);
                    pm1->see(-1, i);
                    pm1->stream(logger, 0);
                    logger << " see " << i << " => ";
                    pm1->stream(logger, -1);
                    logger << std::endl;
                    if (pm1->compare(-1, 1) < 0) {
                        logger << "ERROR!" << std::endl;
                        return;
                    }
                    pm1->copy(-1, 1);
                    pm1->copy(0, -3);
                    pm1->inc(0); // Inc marks as top too quickly
                    if (pm1->compare(0, -3) <= 0) {
                        logger << "ERROR2!" << std::endl;
                        return;
                    }
                }
            }
            return;
        }
    }   
}
