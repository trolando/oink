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

# Vertex IDs differ from their positions in priority order. The error must
# refer to the original IDs, and --sol must reject the file even without -v.
file(WRITE "${WORK}/closure.pg" "parity 4;\n0 2 0 0,2;\n1 4 1 0,2;\n2 1 1 2;\n3 0 0 3;\n")
file(WRITE "${WORK}/closed.sol" "0 0 0;\n")
run_oink("closed partial solution" -v --tl --sol "${WORK}/closed.sol" "${WORK}/closure.pg")

function(reject_solution name solution expected)
    file(WRITE "${WORK}/open.sol" "${solution}")
    execute_process(
        COMMAND "${OINK}" --tl --sol "${WORK}/open.sol" "${WORK}/closure.pg"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
        ERROR_VARIABLE errors
        TIMEOUT 60
    )
    if(result EQUAL 0 OR NOT "${output}${errors}" MATCHES "${expected}")
        message(FATAL_ERROR "${name} gave an unexpected result (${result}):\n${output}${errors}")
    endif()
endfunction()

reject_solution("opponent can escape" "0 0 0;\n1 0;\n"
                "invalid solution: loser can escape on move from vertex 1 to vertex 2")
reject_solution("strategy to an unsolved vertex" "0 0 2;\n"
                "invalid solution: strategy leaves dominion on move from vertex 0 to vertex 2")
reject_solution("strategy to the other region" "0 0 2;\n2 1 2;\n"
                "invalid solution: strategy leaves dominion on move from vertex 0 to vertex 2")
reject_solution("multi-strategy with an escape" "0 0 0,2;\n"
                "invalid solution: strategy leaves dominion on move from vertex 0 to vertex 2")
reject_solution("strategy outside the game" "0 0 9;\n"
                "invalid solution: strategy is not a valid move from vertex 0 to vertex 9")
