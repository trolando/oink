# Copyright 2026 Tom van Dijk
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# The region-renumber test of test_unit on the counter families, on which the Zielonka solvers
# take thousands of steps with a deep stack. TOOLS holds the paths of the generators, separated
# by "|"; the games are written to WORK.
string(REPLACE "|" ";" tools "${TOOLS}")
file(MAKE_DIRECTORY "${WORK}")
set(games "")
foreach(tool ${tools})
    get_filename_component(name "${tool}" NAME_WE)
    set(sizes 4 6 8)
    if(name STREQUAL "counter_rob")
        set(sizes 4 6) # (ZLKW takes 100,000 steps on counter_rob 8, slow with the sanitizers)
    endif()
    foreach(n ${sizes})
        set(game "${WORK}/${name}-${n}.pg")
        execute_process(
            COMMAND "${tool}" ${n}
            OUTPUT_FILE "${game}"
            RESULT_VARIABLE result
        )
        if(NOT result EQUAL 0)
            message(FATAL_ERROR "Could not make ${name} ${n} (${result})")
        endif()
        list(APPEND games "${game}")
    endforeach()
endforeach()
execute_process(
    COMMAND "${TEST_UNIT}" region-renumber-files ${games}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE errors
)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "The region-renumber test failed (${result}):\n${output}${errors}")
endif()
message("${output}")
