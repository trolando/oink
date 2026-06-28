cmake_minimum_required(VERSION 3.14)

project(oink
    VERSION 1.0.0
    DESCRIPTION "Oink, an implementation of modern parity game solvers"
    HOMEPAGE_URL "https://github.com/trolando/oink"
    LANGUAGES C CXX
)

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
        GIT_TAG        v2.3.2
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
    include/oink/player.hpp
    include/oink/solution.hpp
    include/oink/multistrategy.hpp
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
    src/game_builder.cpp
    src/io.cpp
    src/pgparser.cpp
    src/oink.cpp
    src/scc.cpp
    src/solvers.cpp
    src/solver.cpp
    src/verifier.cpp
    # Solvers
    src/solvers/zlk.cpp
    src/solvers/zlkq.cpp
    src/solvers/npp.cpp
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
    src/solvers/rtl.cpp
    src/solvers/zlkpp.cpp
    src/solvers/ptl.cpp
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

    if(OINK_BUILD_TOOLS)
        add_executable(test_solvers)
        target_sources(test_solvers PRIVATE test/test_solvers.cpp)
        set_target_props(test_solvers)
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
    endif()

    if(OINK_BUILD_TESTS)
        enable_testing()
        set(CTEST_TESTING_TIMEOUT 30)

        # parser unit tests (malformed and valid input)
        add_executable(test_parser)
        target_sources(test_parser PRIVATE test/test_parser.cpp)
        set_target_props(test_parser)
        add_test(NAME TestParser COMMAND test_parser)

        # Solution class unit tests
        add_executable(test_solution)
        target_sources(test_solution PRIVATE test/test_solution.cpp)
        set_target_props(test_solution)
        add_test(NAME TestSolution COMMAND test_solution)

        # MultiStrategy class unit tests
        add_executable(test_multistrategy)
        target_sources(test_multistrategy PRIVATE test/test_multistrategy.cpp)
        set_target_props(test_multistrategy)
        add_test(NAME TestMultiStrategy COMMAND test_multistrategy)

        # GameBuilder unit tests
        add_executable(test_builder)
        target_sources(test_builder PRIVATE test/test_builder.cpp)
        set_target_props(test_builder)
        add_test(NAME TestBuilder COMMAND test_builder)

        # Game unit tests (dynamic growth)
        add_executable(test_game)
        target_sources(test_game PRIVATE test/test_game.cpp)
        set_target_props(test_game)
        add_test(NAME TestGame COMMAND test_game)

        # Verifier unit tests (must reject bad strategies)
        add_executable(test_verifier)
        target_sources(test_verifier PRIVATE test/test_verifier.cpp)
        set_target_props(test_verifier)
        add_test(NAME TestVerifier COMMAND test_verifier)

        # test fixpoint algorithms
        add_test(NAME TestSolverFPIpar COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpi -w 0)
        add_test(NAME TestSolverFPIseq COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpi -w -1)
        add_test(NAME TestSolverFPIM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpim)
        add_test(NAME TestSolverFPJ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpj)
        add_test(NAME TestSolverFPJG COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpjg)
        add_test(NAME TestSolverFPJM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpjm)
        # test strategy iteration
        add_test(NAME TestSolverPSIs COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --psi -w -1)
        add_test(NAME TestSolverPSIp COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --psi -w 0)
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
        add_test(NAME TestSolverPPP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ppp)
        add_test(NAME TestSolverRR COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rr)
        add_test(NAME TestSolverDP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --dp)
        add_test(NAME TestSolverRRDP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rrdp)
        add_test(NAME TestSolverPPQ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ppq)
        # test TL and variations
        add_test(NAME TestSolverTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tl)
        add_test(NAME TestSolverRTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rtl)
        add_test(NAME TestSolverORTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ortl)
        add_test(NAME TestSolverPTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ptl)
        add_test(NAME TestSolverSPPTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --spptl)
        add_test(NAME TestSolverDTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --dtl)
        add_test(NAME TestSolverIDTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --idtl)
        add_test(NAME TestSolverTLQ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tlq)
        add_test(NAME TestSolverGPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --gpm)
        add_test(NAME TestSolverPMTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pmtl)
        # test ZLK variations
        add_test(NAME TestSolverZLKseq COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlk -w -1)
        add_test(NAME TestSolverZLKpar COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlk -w 0)
        add_test(NAME TestSolverZLKQ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkq)
        add_test(NAME TestSolverUZLK COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --uzlk)
        add_test(NAME TestSolverZLKPP1 COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkpp-std)
        add_test(NAME TestSolverZLKPP2 COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkpp-waw)
        add_test(NAME TestSolverZLKPP3 COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkpp-liv)
        # Template for adding a new solver test (replace -- with the solver flag):
        # add_test(NAME TestSolver COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --)
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

