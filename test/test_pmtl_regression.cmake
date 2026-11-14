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

# The sparse seeded random game on which PMTL with small measures counted
# through the cycles of known tangles instead of lifting their heads to the
# escape measures. Only test_solvers is used, which also builds on Windows.
execute_process(
    COMMAND "${TEST_SOLVERS}" --gameseed=4066755396 --size=1500
    OUTPUT_FILE "${GAME}"
    RESULT_VARIABLE result
)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Could not generate the PMTL-small regression game (${result})")
endif()
execute_process(
    COMMAND "${TEST_SOLVERS}" --pmtl-small "${GAME}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE errors
    TIMEOUT 20
)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "PMTL-small regression failed (${result}):\n${output}${errors}")
endif()
