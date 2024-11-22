cmake_minimum_required(VERSION 3.14)

project(oink
    VERSION 1.0.0
    DESCRIPTION "Oink, an implementation of modern parity game solvers"
    HOMEPAGE_URL "https://github.com/trolando/oink"
    LANGUAGES C CXX
)

# Dependencies
include(FetchContent)
if(NOT TARGET lace)
  find_package(lace QUIET)
  if(NOT lace_FOUND)
    FetchContent_Declare(
        lace
        GIT_REPOSITORY https://github.com/trolando/lace.git
        GIT_TAG        v1.4.1
    )
    FetchContent_MakeAvailable(lace)
  endif()
endif()

find_package(Boost REQUIRED COMPONENTS filesystem iostreams regex system random CONFIG)

# Main library target
add_library(oink STATIC)
add_library(oink::oink ALIAS oink)
set_target_properties(oink PROPERTIES VERSION ${oink_VERSION} SOVERSION ${oink_VERSION_MAJOR})

set(OINK_HDRS
    include/oink/oink.hpp
    include/oink/error.hpp
    include/oink/game.hpp
    include/oink/solvers.hpp
    include/oink/solver.hpp
    include/oink/bitset.hpp
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
    $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(oink PUBLIC pthread lace::lace)

option(OINK_BUILD_NATIVE "Use -march=native flag for compiling Oink" ON)
if(OINK_BUILD_NATIVE)
    target_compile_options(oink PRIVATE "-march=native")
endif()

target_sources(oink
  PRIVATE
    # Core files
    src/game.cpp
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
    src/solvers/fpj.cpp
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
    ${OINK_HDRS}
)

set_target_properties(oink PROPERTIES PUBLIC_HEADER "${OINK_HDRS}")

# Common target properties
macro(set_target_props NAME)
    target_compile_features(${NAME} PUBLIC c_std_17 cxx_std_17)
    target_compile_options(${NAME} PRIVATE -Wall -Wextra)
    target_include_directories(${NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${Boost_INCLUDE_DIRS})
    target_link_libraries(${NAME} oink::oink)
endmacro(set_target_props)

add_executable(solve)
target_sources(solve PRIVATE src/tools/solve.cpp src/tools/getrss.c)
set_target_properties(solve PROPERTIES OUTPUT_NAME oink)
set_target_props(solve)
target_link_libraries(solve Boost::iostreams Boost::random)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(OINK_BUILD_TOOLS "Build basic Oink tools" ON)
    option(OINK_BUILD_EXTRA_TOOLS "Build extra Oink tools" OFF)
    option(OINK_BUILD_TESTS "Build tests" ON)

    if(OINK_BUILD_TOOLS)
        add_executable(test_solvers)
        target_sources(test_solvers PRIVATE test/test_solvers.cpp)
        set_target_props(test_solvers)
        target_link_libraries(test_solvers Boost::iostreams Boost::filesystem Boost::random)
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

        # test fixpoint algorithms
        add_test(NAME TestSolverFPIpar COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpi -w 0)
        add_test(NAME TestSolverFPIseq COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpi -w -1)
        add_test(NAME TestSolverFPJ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpj)
        add_test(NAME TestSolverFPJG COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpjg)
        # test strategy iteration
        add_test(NAME TestSolverPSIs COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --psi -w -1)
        add_test(NAME TestSolverPSIp COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --psi -w 0)
        add_test(NAME TestSolverSSI COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ssi)
        # test progress measures and variations
        add_test(NAME TestSolverSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --spm)
        add_test(NAME TestSolverTSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tspm)
        # THERE IS A BUG IN MSPM -- disable test for now, but it needs to be solved
        #  add_test(NAME TestSolverMSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --mspm)
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
        # THERE IS A BUG IN RRDP -- disable test for now, but it needs to be solved
        # add_test(NAME TestSolverRRDP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rrdp)
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
        # test ZLK variations
        add_test(NAME TestSolverZLKseq COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlk -w -1)
        add_test(NAME TestSolverZLKpar COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlk -w 0)
        add_test(NAME TestSolverZLKQ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkq)
        add_test(NAME TestSolverUZLK COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --uzlk)
        add_test(NAME TestSolverZLKPP1 COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkpp-std)
        add_test(NAME TestSolverZLKPP2 COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkpp-waw)
        add_test(NAME TestSolverZLKPP3 COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkpp-liv)
        #add_test(NAME TestSolver COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --)
    endif()
endif()

set_target_properties(oink PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
)

include(GNUInstallDirs)

install(TARGETS oink solve
    EXPORT oink-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oink
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/oink
)

install(EXPORT oink-targets
    FILE oink-targets.cmake
    NAMESPACE oink::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/oink
)

include(CMakePackageConfigHelpers)

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

export(EXPORT oink-targets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/oink-targets.cmake
    NAMESPACE oink::
)
