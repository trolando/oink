cmake_minimum_required(VERSION 3.14)

project(oink
    VERSION 1.2.1
    DESCRIPTION "Oink, an implementation of modern parity game solvers"
    HOMEPAGE_URL "https://github.com/trolando/oink"
    LANGUAGES C CXX
)

# Default to an optimized build when no build type is given. Only for
# single-config generators, and only when Oink is the top-level project, so a
# parent project's choice is never overridden.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR AND NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type given, defaulting to Release")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
endif()

# Build Oink (and the bundled lace) tuned to the build machine's CPU. Off by
# default so packaged binaries stay portable; enable for local benchmarking.
option(OINK_BUILD_NATIVE "Optimize for the build CPU (-march=native; not portable)" OFF)

# Dependencies
include(FetchContent)
if(NOT TARGET lace::lace)
  find_package(lace QUIET)
  if(NOT lace_FOUND)
    FetchContent_Declare(
        lace
        GIT_REPOSITORY https://github.com/trolando/lace.git
        GIT_TAG        v2.6.0
    )
    # Keep the bundled lace's CPU tuning consistent with Oink's. Uses the cache
    # default so an explicit -DLACE_NATIVE_OPT still wins.
    set(LACE_NATIVE_OPT ${OINK_BUILD_NATIVE} CACHE BOOL "Enable -march=native in lace")
    FetchContent_MakeAvailable(lace)
  endif()
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Main library target
add_library(oink STATIC)
add_library(oink::oink ALIAS oink)

# Generate the version header from the single source of truth: project(VERSION).
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/oink/version.hpp
    @ONLY
)

set(OINK_HDRS
    include/oink/oink.hpp
    include/oink/error.hpp
    include/oink/game.hpp
    include/oink/game_builder.hpp
    include/oink/io.hpp
    include/oink/pgparser.hpp
    include/oink/player.hpp
    include/oink/solution.hpp
    include/oink/span.hpp
    include/oink/subgame_view.hpp
    include/oink/solvers.hpp
    include/oink/solver.hpp
    include/oink/bitset.hpp
    include/oink/intrinsics.hpp
    include/oink/uintqueue.hpp
    include/oink/libpopcnt.h
)

target_compile_features(oink PUBLIC c_std_17 cxx_std_17)

target_compile_options(oink PRIVATE
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:
        -pipe
        $<$<CONFIG:Debug>:-O0 -Wall -Wextra>>
    $<$<CXX_COMPILER_ID:MSVC>:
        $<$<CONFIG:Debug>:/Od /Wall /Zi>>
)

target_include_directories(oink
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/src
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(oink PUBLIC Threads::Threads lace::lace)

# Optional transparent decompression of input parity games (.gz/.bz2/.xz).
# These are private implementation details of the parser: each backend is
# independent, a missing library simply disables that format (reported as an
# informational message at configure time), and the dependency does not appear
# in Oink's public headers.
# OINK_COMPRESSION_DEPS feeds find_dependency() in the CMake package config;
# OINK_PC_EXTRA_LIBS feeds Libs.private in the pkg-config file.
set(OINK_COMPRESSION_DEPS "")
set(OINK_PC_EXTRA_LIBS "")
find_package(ZLIB)
if(ZLIB_FOUND)
    target_link_libraries(oink PRIVATE ZLIB::ZLIB)
    target_compile_definitions(oink PRIVATE OINK_HAVE_ZLIB=1)
    list(APPEND OINK_COMPRESSION_DEPS ZLIB)
    string(APPEND OINK_PC_EXTRA_LIBS " -lz")
    message(STATUS "Oink: gzip (.gz) input support enabled")
else()
    message(STATUS "Oink: gzip (.gz) input support disabled (zlib not found)")
endif()

find_package(BZip2)
if(BZIP2_FOUND)
    target_link_libraries(oink PRIVATE BZip2::BZip2)
    target_compile_definitions(oink PRIVATE OINK_HAVE_BZIP2=1)
    list(APPEND OINK_COMPRESSION_DEPS BZip2)
    string(APPEND OINK_PC_EXTRA_LIBS " -lbz2")
    message(STATUS "Oink: bzip2 (.bz2) input support enabled")
else()
    message(STATUS "Oink: bzip2 (.bz2) input support disabled (libbz2 not found)")
endif()

find_package(LibLZMA)
if(LibLZMA_FOUND)
    target_link_libraries(oink PRIVATE LibLZMA::LibLZMA)
    target_compile_definitions(oink PRIVATE OINK_HAVE_LZMA=1)
    list(APPEND OINK_COMPRESSION_DEPS LibLZMA)
    string(APPEND OINK_PC_EXTRA_LIBS " -llzma")
    message(STATUS "Oink: xz (.xz) input support enabled")
else()
    message(STATUS "Oink: xz (.xz) input support disabled (liblzma not found)")
endif()

# OINK_BUILD_NATIVE is declared near the top so it can be propagated to lace.
if(OINK_BUILD_NATIVE)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-march=native" OINK_HAS_MARCH_NATIVE)
    if(OINK_HAS_MARCH_NATIVE)
        target_compile_options(oink PRIVATE -march=native)
    else()
        message(WARNING "OINK_BUILD_NATIVE is ON but the compiler does not support -march=native")
    endif()
endif()

option(OINK_ENABLE_SANITIZERS "Enable address and undefined behavior sanitizers" OFF)
if(OINK_ENABLE_SANITIZERS)
    target_compile_options(oink PUBLIC
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
            -g>
    )
    target_link_options(oink PUBLIC
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
            -fsanitize=address,undefined>
    )
endif()

# ThreadSanitizer is a separate option: it detects data races (which ASan/UBSan do
# not) and cannot be combined with them. Use this to verify the parallel solvers.
option(OINK_ENABLE_TSAN "Enable the thread sanitizer (data-race detection)" OFF)
if(OINK_ENABLE_TSAN)
    if(OINK_ENABLE_SANITIZERS)
        message(FATAL_ERROR "OINK_ENABLE_TSAN cannot be combined with OINK_ENABLE_SANITIZERS")
    endif()
    target_compile_options(oink PUBLIC
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
            -fsanitize=thread
            -fno-omit-frame-pointer
            -g>
    )
    target_link_options(oink PUBLIC
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
            -fsanitize=thread>
    )
endif()

option(OINK_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)
if(OINK_WARNINGS_AS_ERRORS)
    target_compile_options(oink PRIVATE
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:-Werror>
        $<$<CXX_COMPILER_ID:MSVC>:/WX>
    )
endif()

target_sources(oink
  PRIVATE
    # Core files
    src/game.cpp
    src/solution.cpp
    src/game_builder.cpp
    src/bzip2_parallel.cpp
    src/io.cpp
    src/pgparser.cpp
    src/oink.cpp
    src/scc.cpp
    src/solvers.cpp
    src/solver.cpp
    src/verifier.cpp
    # Solvers
    src/solvers/zlk.cpp
    src/solvers/zlkj.cpp
    src/solvers/zlkq.cpp
    src/solvers/npp.cpp
    src/solvers/ppj.cpp
    src/solvers/pp.cpp
    src/solvers/ppp.cpp
    src/solvers/rr.cpp
    src/solvers/dp.cpp
    src/solvers/rrdp.cpp
    src/solvers/ppq.cpp
    src/solvers/fpi.cpp
    src/solvers/fpim.cpp
    src/solvers/fpj.cpp
    src/solvers/fpjm.cpp
    src/solvers/psi.cpp
    src/solvers/ssi.cpp
    src/solvers/spm.cpp
    src/solvers/mspm.cpp
    src/solvers/tspm.cpp
    src/solvers/sspm.cpp
    src/solvers/qpt.cpp
    src/solvers/tl.cpp
    src/solvers/tlpar.cpp
    src/solvers/tlj.cpp
    src/solvers/rtl.cpp
    src/solvers/zlkpp.cpp
    src/solvers/dtl.cpp
    src/solvers/tlq.cpp
    src/solvers/gpm.cpp
    src/solvers/pmtl.cpp
    ${OINK_HDRS}
)

set_target_properties(oink PROPERTIES PUBLIC_HEADER "${OINK_HDRS}")

# Common target properties
macro(set_target_props NAME)
    target_compile_features(${NAME} PUBLIC c_std_17 cxx_std_17)
    target_compile_options(${NAME} PRIVATE -Wall -Wextra)
    target_include_directories(${NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(${NAME} oink::oink)
endmacro(set_target_props)

# Command-line tools. On by default; downstream consumers that embed Oink as a
# library (e.g. via add_subdirectory) can disable them with -DOINK_BUILD_TOOLS=OFF.
option(OINK_BUILD_TOOLS "Build the Oink command-line tools" ON)

if(OINK_BUILD_TOOLS)
    add_executable(solve)
    target_sources(solve PRIVATE src/tools/solve.cpp src/tools/getrss.c)
    set_target_properties(solve PROPERTIES OUTPUT_NAME oink)
    set_target_props(solve)
endif()

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(OINK_BUILD_EXTRA_TOOLS "Build extra Oink tools" OFF)
    option(OINK_BUILD_EXAMPLES "Build the library usage examples" ON)
    option(OINK_BUILD_TESTS "Build tests" ON)

    # Support code shared by the test tools and their unit tests: the terminal,
    # the interactive solver chooser and external solvers.
    if(OINK_BUILD_TOOLS OR OINK_BUILD_TESTS)
        add_library(oink_testtools STATIC
            test/common/terminal.cpp
            test/common/chooser.cpp
            test/common/external_solver.cpp
        )
        set_target_props(oink_testtools)

        # bench_solvers, apart from its main, so the unit tests can use it
        add_library(oink_bench STATIC
            test/bench/core.cpp
            test/bench/display.cpp
            test/bench/experiment.cpp
            test/bench/live.cpp
            test/bench/options.cpp
            test/bench/report.cpp
            test/bench/runner.cpp
            test/bench/setup.cpp
            test/bench/summary.cpp
            test/bench/table.cpp
            test/bench/trial.cpp
            test/bench/worker.cpp
        )
        target_include_directories(oink_bench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test)
        target_compile_definitions(oink_bench PRIVATE OINK_BENCH_BUILD_TYPE="$<CONFIG>")
        target_link_libraries(oink_bench oink_testtools)
        set_target_props(oink_bench)
    endif()

    if(OINK_BUILD_TOOLS)
        add_executable(test_solvers)
        target_sources(test_solvers PRIVATE test/test_solvers.cpp)
        target_link_libraries(test_solvers oink_testtools)
        set_target_props(test_solvers)
        add_executable(bench_solvers test/bench/main.cpp)
        target_include_directories(bench_solvers PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test)
        target_link_libraries(bench_solvers oink_bench oink_testtools)
        set_target_props(bench_solvers)
        add_executable(test_measures)
        target_sources(test_measures PRIVATE test/test_measures.cpp)
        set_target_props(test_measures)
    endif()

    # A small, self-contained example of using Oink as a library. Built by
    # default so it keeps compiling against the public API; see examples/simple.cpp.
    if(OINK_BUILD_EXAMPLES)
        add_executable(oink-example-simple examples/simple.cpp)
        set_target_props(oink-example-simple)
    endif()

    if(OINK_BUILD_EXTRA_TOOLS)
        add_executable(nudge src/tools/nudge.cpp)
        set_target_props(nudge)

        add_executable(dotty src/tools/dotty.cpp)
        set_target_props(dotty)

        add_executable(verify src/tools/verify.cpp)
        set_target_props(verify)

        add_executable(simple src/tools/simple.cpp)
        set_target_props(simple)

        add_executable(rngame src/tools/rngame.cpp)
        set_target_props(rngame)

        add_executable(stgame src/tools/stgame.cpp)
        set_target_props(stgame)

        add_executable(counter_rr src/tools/counter_rr.cpp)
        set_target_props(counter_rr)

        add_executable(counter_dp src/tools/counter_dp.cpp)
        set_target_props(counter_dp)

        add_executable(counter_m src/tools/counter_m.cpp)
        set_target_props(counter_m)

        add_executable(counter_core src/tools/counter_core.cpp)
        set_target_props(counter_core)

        add_executable(counter_rob src/tools/counter_rob.cpp)
        set_target_props(counter_rob)

        add_executable(counter_symsi src/tools/counter_symsi.cpp)
        set_target_props(counter_symsi)

        add_executable(counter_ortl src/tools/counter_ortl.cpp)
        set_target_props(counter_ortl)

        add_executable(counter_qpt src/tools/counter_qpt.cpp)
        set_target_props(counter_qpt)

        add_executable(tc src/tools/tc.cpp)
        set_target_props(tc)

        add_executable(tc+ src/tools/tc+.cpp)
        set_target_props(tc+)

        add_executable(multistrat src/tools/multistrat.cpp)
        set_target_props(multistrat)
    endif()

    if(OINK_BUILD_TESTS)
        enable_testing()
        set(CTEST_TESTING_TIMEOUT 30)

        # One unit-test runner; CTest still reports each check separately.
        add_executable(test_unit
            test/test_main.cpp
            test/test_measure_ops.cpp
            test/test_parser.cpp
            test/test_io.cpp
            test/test_solution.cpp
            test/test_ppj.cpp
            test/test_zlkj.cpp
            test/test_fpj_batch.cpp
            test/test_pma.cpp
            test/test_fpi_frontier.cpp
            test/test_tlq_cache.cpp
            test/test_tlj_cache.cpp
            test/test_gpm_cache.cpp
            test/test_pmtl_cache.cpp
            test/test_in_edge_map.cpp
            test/test_builder.cpp
            test/test_game.cpp
            test/test_verifier.cpp
            test/test_terminal.cpp
            test/test_chooser.cpp
            test/test_bench.cpp
            test/test_table.cpp
        )
        target_include_directories(test_unit PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test)
        target_link_libraries(test_unit oink_bench oink_testtools)
        set_target_props(test_unit)
        target_compile_definitions(test_unit PRIVATE OINK_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests")
        if(BZIP2_FOUND)
            target_link_libraries(test_unit BZip2::BZip2)
            target_compile_definitions(test_unit PRIVATE OINK_HAVE_BZIP2=1)
        endif()

        add_test(NAME TestMeasureOps COMMAND test_unit measure-ops)
        add_test(NAME TestParser COMMAND test_unit parser)
        set(parser_input "${CMAKE_CURRENT_SOURCE_DIR}/test/parser-input/records.pg")
        foreach(parser fast renumber)
            add_test(NAME TestParserStdin-${parser}
                COMMAND ${CMAKE_COMMAND} "-DTEST_UNIT=$<TARGET_FILE:test_unit>"
                    "-DMODE=parser-input-${parser}" "-DINPUT=${parser_input}"
                    -P "${CMAKE_CURRENT_SOURCE_DIR}/test/test_parser_stdin.cmake")
            if(ZLIB_FOUND)
                add_test(NAME TestParserGzip-${parser}
                    COMMAND test_unit parser-input-${parser} "${parser_input}.gz")
            endif()
            if(BZIP2_FOUND)
                add_test(NAME TestParserBzip2-${parser}
                    COMMAND test_unit parser-input-${parser} "${parser_input}.bz2")
            endif()
            if(LibLZMA_FOUND)
                add_test(NAME TestParserXz-${parser}
                    COMMAND test_unit parser-input-${parser} "${parser_input}.xz")
            endif()
        endforeach()
        add_test(NAME TestIO COMMAND test_unit io)
        add_test(NAME TestSolution COMMAND test_unit solution)
        add_test(NAME TestFPJBatch COMMAND test_unit fpj-batch)
        add_test(NAME TestPPJ COMMAND test_unit ppj)
        add_test(NAME TestZLKJ COMMAND test_unit zlkj)
        add_test(NAME TestPMA COMMAND test_unit pma)
        add_test(NAME TestFPIFrontier COMMAND test_unit fpi-frontier)
        add_test(NAME TestFPIParallel COMMAND test_unit fpi-parallel)
        add_test(NAME TestTLJCache COMMAND test_unit tlj-cache)
        add_test(NAME TestTLQCache COMMAND test_unit tlq-cache)
        add_test(NAME TestGPMCache COMMAND test_unit gpm-cache)
        add_test(NAME TestPMTLCache COMMAND test_unit pmtl-cache)
        add_test(NAME TestInEdgeMap COMMAND test_unit in-edge-map)
        add_test(NAME TestBuilder COMMAND test_unit builder)
        add_test(NAME TestGame COMMAND test_unit game)
        add_test(NAME TestVerifier COMMAND test_unit verifier)
        add_test(NAME TestTerminal COMMAND test_unit terminal)
        add_test(NAME TestChooser COMMAND test_unit chooser)
        add_test(NAME TestBench COMMAND test_unit bench)
        add_test(NAME TestTable COMMAND test_unit table)

        if(OINK_BUILD_TOOLS) # the solver tests use test_solvers
            # test fixpoint algorithms
            add_test(NAME TestSolverFPIpar COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpi -w 0 --force-parallel)
            add_test(NAME TestSolverFPIseq COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpi -w -1)
            add_test(NAME TestSolverFPIM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpim)
            add_test(NAME TestSolverFPJ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpj)
            add_test(NAME TestSolverFPJG COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpjg)
            add_test(NAME TestSolverFPJM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpjm)
            # test strategy iteration
            add_test(NAME TestSolverPSIs COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --psi -w -1)
            add_test(NAME TestSolverPSIp COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --psi -w 0 --force-parallel -c parmin=1,chunk=4)
            add_test(NAME TestSolverSSI COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ssi)
            # test progress measures and variations
            add_test(NAME TestSolverSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --spm)
            add_test(NAME TestSolverTSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tspm)
            add_test(NAME TestSolverMSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --mspm)
            add_test(NAME TestSolverSSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --sspm)
            add_test(NAME TestSolverQPT COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --qpt)
            add_test(NAME TestSolverBSSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --bsspm)
            add_test(NAME TestSolverBQPT COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --bqpt)
            # test PP and variations
            add_test(NAME TestSolverNPP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --npp)
            add_test(NAME TestSolverPP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pp)
            add_test(NAME TestSolverPPJ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ppj)
            add_test(NAME TestSolverPPJscc COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ppj --scc)
            add_test(NAME TestSolverPPP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ppp)
            add_test(NAME TestSolverRR COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rr)
            add_test(NAME TestSolverDP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --dp)
            add_test(NAME TestSolverRRDP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rrdp)
            add_test(NAME TestSolverPPQ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ppq)
            # test TL and variations
            add_test(NAME TestSolverTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tl)
            add_test(NAME TestSolverTLPARseq COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tlpar -w -1)
            add_test(NAME TestSolverTLPAR1 COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tlpar -w 1 --force-parallel)
            add_test(NAME TestSolverTLPARpar COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tlpar -w 0 --force-parallel -c parmin=1,chunk=4)
            add_test(NAME TestSolverTLPARcfg COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tlpar --tlpar=parmin=1,chunk=4 -w 0 --force-parallel)
            add_test(NAME TestSolverTLJ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tlj)
            add_test(NAME TestSolverRTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rtl)
            add_test(NAME TestSolverORTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ortl)
            add_test(NAME TestSolverDTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --dtl)
            add_test(NAME TestSolverIDTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --idtl)
            add_test(NAME TestSolverTLQ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tlq)
            add_test(NAME TestSolverGPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --gpm)
            add_test(NAME TestSolverPMTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pmtl)
            add_test(NAME TestSolverGPMsmall COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --gpm-small)
            add_test(NAME TestSolverGPMsuccinct COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --gpm-succinct)
            add_test(NAME TestSolverGPMordered COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --gpm-ordered)
            add_test(NAME TestSolverPMTLsmall COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pmtl-small)
            add_test(NAME TestSolverPMTLsmallRegression
                COMMAND ${CMAKE_COMMAND}
                    -DTEST_SOLVERS=$<TARGET_FILE:test_solvers>
                    -DGAME=${CMAKE_CURRENT_BINARY_DIR}/pmtl-small-regression.pg
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test/test_pmtl_regression.cmake)
            add_test(NAME TestSolverPMTLsuccinct COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pmtl-succinct)
            add_test(NAME TestSolverPMTLordered COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pmtl-ordered)
            add_test(NAME TestSolverGPMstrahler COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --gpm-strahler)
            add_test(NAME TestSolverPMTLstrahler COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pmtl-strahler)
            add_test(NAME TestSolverPMA COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pma)
            add_test(NAME TestSolverPMAsmall COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pma-small)
            add_test(NAME TestSolverPMAsuccinct COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pma-succinct)
            add_test(NAME TestSolverPMAordered COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pma-ordered)
            add_test(NAME TestSolverPMAstrahler COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pma-strahler)
            add_test(NAME TestSolverPMAscc COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pma --scc)
            # test ZLK variations
            add_test(NAME TestSolverZLKseq COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlk -w -1)
            add_test(NAME TestSolverZLKpar COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlk -w 0 --force-parallel)
            add_test(NAME TestSolverZLKQ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkq)
            add_test(NAME TestSolverUZLK COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --uzlk)
            add_test(NAME TestSolverZLKJ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkj)
            add_test(NAME TestSolverZLKPP1 COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkpp-std)
            add_test(NAME TestSolverZLKPP2 COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkpp-waw)
            add_test(NAME TestSolverZLKPP3 COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkpp-liv)
            # test solving one bottom SCC at a time
            add_test(NAME TestSolverFPIscc COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpi --scc)
            add_test(NAME TestSolverTLscc COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tl --scc)
            add_test(NAME TestSolverZLKscc COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlk --scc)
            add_test(NAME TestSolverZLKJscc COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkj --scc)
            add_test(NAME TestSolverPPscc COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pp --scc)
            # multi-strategies in solution files, with oink and --sol (see test/test_multisol.cmake)
            add_test(NAME TestMultiStrategyFiles
                COMMAND ${CMAKE_COMMAND}
                    -DOINK=$<TARGET_FILE:solve>
                    -DGAME=${CMAKE_CURRENT_SOURCE_DIR}/examples/amba_decomposed_arbiter_4.tlsf.ehoa.pg
                    -DWORK=${CMAKE_CURRENT_BINARY_DIR}/multisol
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test/test_multisol.cmake)
            # short headless runs of bench_solvers (see test/bench/smoke.cmake)
            set(BENCH_MODES basic sampling reference random errors preprocess live experiment)
            if(NOT WIN32)
                list(APPEND BENCH_MODES timeout fresh skip) # worker processes
                # the interactive screens on a pseudo terminal
                add_executable(test_pty test/bench/test_pty.cpp)
                target_include_directories(test_pty PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test)
                target_link_libraries(test_pty oink_testtools)
                set_target_props(test_pty)
                add_test(NAME TestBenchTerminal COMMAND test_pty $<TARGET_FILE:bench_solvers> $<TARGET_FILE:solve> ${CMAKE_CURRENT_SOURCE_DIR}/tests ${CMAKE_CURRENT_BINARY_DIR}/bench-smoke/terminal)
                # several sessions of a few seconds each; slow under the sanitizers
                set_tests_properties(TestBenchTerminal PROPERTIES TIMEOUT 300)
            endif()
            foreach(mode ${BENCH_MODES})
                add_test(NAME TestBench_${mode}
                    COMMAND ${CMAKE_COMMAND}
                        -DBENCH=$<TARGET_FILE:bench_solvers>
                        -DOINK=$<TARGET_FILE:solve>
                        -DDATA=${CMAKE_CURRENT_SOURCE_DIR}/tests
                        -DWORK=${CMAKE_CURRENT_BINARY_DIR}/bench-smoke/${mode}
                        -DMODE=${mode}
                        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/bench/smoke.cmake)
                set_tests_properties(TestBench_${mode} PROPERTIES TIMEOUT 300)
            endforeach()
            # without a terminal, test_solvers -i stops with an error instead of waiting
            add_test(NAME TestSolversInteractiveNoTerminal COMMAND test_solvers -i --tl ${CMAKE_CURRENT_SOURCE_DIR}/tests/vb001)
            set_tests_properties(TestSolversInteractiveNoTerminal PROPERTIES PASS_REGULAR_EXPRESSION "needs a (POSIX )?terminal")
            # Template for adding a new solver test (replace -- with the solver flag):
            # add_test(NAME TestSolver COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --)
        endif()
    endif()
endif()

set_target_properties(oink PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
)

include(GNUInstallDirs)

install(TARGETS oink
    EXPORT oink-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oink
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oink
)

# The generated version header is not part of PUBLIC_HEADER; install it directly.
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/generated/oink/version.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oink
)

if(OINK_BUILD_TOOLS)
    install(TARGETS solve
        EXPORT oink-targets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()

install(EXPORT oink-targets
    FILE oink-targets.cmake
    NAMESPACE oink::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oink
)

include(CMakePackageConfigHelpers)

# oink-config.cmake.in resolves Threads and lace unconditionally and loops over
# OINK_COMPRESSION_DEPS (set during backend detection above) for the optional
# compression libraries that were enabled in this build.
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/oink-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/oink-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oink
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/oink-config-version.cmake
    VERSION ${oink_VERSION}
    COMPATIBILITY SameMinorVersion
)

install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/oink-config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/cmake/oink-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oink
)

# pkg-config file for non-CMake consumers. OINK_PC_EXTRA_LIBS (the enabled
# compression libs) goes into Libs.private for static linking; lace is pulled
# in via Requires.
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/oink.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/oink.pc
    @ONLY
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/oink.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)

