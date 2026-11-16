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

# Multi-strategies in solution files, end to end with oink (OINK) on a game
# whose vertices are not in the order of their priorities (GAME): fpim writes
# its multi-strategy, --sol reads it back, and a partial solution (the region
# of Even) keeps all its moves, with and without solving one SCC at a time.

file(MAKE_DIRECTORY "${WORK}")

function(run_oink name)
    execute_process(
        COMMAND "${OINK}" ${ARGN}
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE errors
        TIMEOUT 60
    )
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "${name} failed (${result}):\n${output}${errors}")
    endif()
endfunction()

run_oink("fpim" --fpim "${GAME}" "${WORK}/fpim.sol")
file(READ "${WORK}/fpim.sol" fpim)
if(NOT fpim MATCHES "[0-9],[0-9]")
    message(FATAL_ERROR "the solution of fpim has no vertex with several moves")
endif()

# a full solution passes through unchanged
run_oink("full solution" -v --tl --sol "${WORK}/fpim.sol" "${GAME}" "${WORK}/again.sol")
file(READ "${WORK}/again.sol" again)
if(NOT again STREQUAL fpim)
    message(FATAL_ERROR "the full solution changed when read and written again")
endif()

# the region of Even as a partial solution (without the semicolons, which
# are optional and separate the items of a list in CMake)
string(REPLACE ";" "" plain "${fpim}")
string(REGEX MATCHALL "\n[0-9]+ 0[^\n]*" lines "\n${plain}")
list(LENGTH lines count)
string(JOIN "" partial ${lines})
file(WRITE "${WORK}/even.sol" "${partial}\n")
foreach(config "--tl" "--tl;--scc" "--fpim" "--fpim;--scc")
    run_oink("partial solution with ${config}" -v ${config} --sol "${WORK}/even.sol" "${GAME}" "${WORK}/out.sol")
    file(READ "${WORK}/out.sol" out)
    string(REPLACE ";" "" out "${out}")
    foreach(line ${lines})
        string(FIND "${out}" "${line}\n" at)
        if(at EQUAL -1)
            message(FATAL_ERROR "with ${config}, the moves of the partial solution changed:${line}")
        endif()
    endforeach()
endforeach()
message(STATUS "${count} vertices of Even kept their moves")
