# SPDX-License-Identifier: Apache-2.0
execute_process(COMMAND "${GENERATOR}" 5 OUTPUT_FILE "${GAME}" RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
    message(FATAL_ERROR "tc+ 5 failed (${result})")
endif()
execute_process(
    COMMAND "${TEST_SOLVERS}" --ep-clipped --ep-sign --gpm-energy --gpa-energy --pmtl-energy "${GAME}"
    RESULT_VARIABLE result TIMEOUT 10
)
if(NOT result STREQUAL "0")
    message(FATAL_ERROR "EP or an energy variant failed on tc+ 5 (${result})")
endif()
