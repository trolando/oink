execute_process(COMMAND "${TEST_UNIT}" "${MODE}" - INPUT_FILE "${INPUT}"
    RESULT_VARIABLE result)
if(NOT result STREQUAL "0")
    message(FATAL_ERROR "Parsing the game from stdin failed (${result})")
endif()
