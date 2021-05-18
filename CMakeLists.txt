cmake_minimum_required(VERSION 3.0)
project(oink VERSION 1.0.0 LANGUAGES C CXX)

set(PROJECT_DESCRIPTION "Oink, an implementation of modern parity game solvers")
set(PROJECT_URL "https://github.com/trolando/oink")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "-Wall -Wextra")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Some info
message(STATUS "OS: ${CMAKE_SYSTEM_NAME}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 4.9
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
        message(FATAL_ERROR "GCC version must be at least 4.9!  " ${CMAKE_CXX_COMPILER_VERSION})
    endif()
endif()

include_directories(src)
include(GNUInstallDirs)

option(BUILD_SHARED_LIBS "Enable/disable creation of shared libraries" OFF)
option(BUILD_STATIC_LIBS "Enable/disable creation of static libraries" ON)
option(BUILD_EXTRA_TOOLS "Enable/disable building of extra tools" OFF)

find_package(Boost REQUIRED filesystem iostreams regex system)

add_library(oink
    # Core files
    src/game.cpp
    src/lace.c
    src/oink.cpp
    src/scc.cpp
    src/solvers.cpp
    src/verifier.cpp
    # Solvers
    src/zlk.cpp
    src/zlkq.cpp
    src/npp.cpp
    src/pp.cpp
    src/ppp.cpp
    src/rr.cpp
    src/dp.cpp
    src/rrdp.cpp
    src/fpi.cpp
    src/fpj.cpp
    src/psi.cpp
    src/spm.cpp
    src/mspm.cpp
    src/tspm.cpp
    src/sspm.cpp
    src/qpt.cpp
    src/tl.cpp
    src/rtl.cpp
    src/zlkpp.cpp
    src/ptl.cpp
    src/dtl.cpp
)
target_link_libraries(oink pthread)

add_executable(solve src/tools/solve.cpp src/tools/getrss.c)
set_target_properties(solve PROPERTIES OUTPUT_NAME oink)
target_link_libraries(solve oink ${Boost_LIBRARIES})

add_executable(test_solvers test/test_solvers.cpp)
target_link_libraries(test_solvers oink ${Boost_LIBRARIES})

enable_testing()

add_test(NAME TestSolverZLKs COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlk -w -1)
add_test(NAME TestSolverZLKp COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlk -w 0)
add_test(NAME TestSolverPP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --pp)
add_test(NAME TestSolverPPP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ppp)
add_test(NAME TestSolverRR COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rr)
add_test(NAME TestSolverDP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --dp)
add_test(NAME TestSolverRRDP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rrdp)
add_test(NAME TestSolverPSIs COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --psi -w -1)
add_test(NAME TestSolverPSIp COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --psi -w 0)
add_test(NAME TestSolverTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tl)
add_test(NAME TestSolverATL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --atl)
add_test(NAME TestSolverOTFTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --otftl)
add_test(NAME TestSolverOTFATL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --otfatl)
add_test(NAME TestSolverRTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rtl)
# add_test(NAME TestSolverORTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ortl)
add_test(NAME TestSolverPTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ptl)
add_test(NAME TestSolverSPPTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --spptl)
add_test(NAME TestSolverDTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --dtl)
add_test(NAME TestSolverIDTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --idtl)
add_test(NAME TestSolverSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --spm)
add_test(NAME TestSolverZLKQs COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkq -w -1)
add_test(NAME TestSolverZLKQp COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --zlkq -w 0)
add_test(NAME TestSolverFPI COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpi)
add_test(NAME TestSolverFPJ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpj)
add_test(NAME TestSolverFPJG COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpjg)
add_test(NAME TestSolverTSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tspm)
#add_test(NAME TestSolverMSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --mspm)
add_test(NAME TestSolverSSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --sspm)
add_test(NAME TestSolverQPT COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --qpt)
add_test(NAME TestSolverBSSPM COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --bsspm)
add_test(NAME TestSolverBQPT COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --bqpt)
if(BUILD_EXTRA_TOOLS)
    add_executable(nudge src/tools/nudge.cpp)
    target_link_libraries(nudge oink)

    add_executable(dotty src/tools/dotty.cpp)
    target_link_libraries(dotty oink)

    add_executable(verify src/tools/verify.cpp)
    target_link_libraries(verify oink)

    add_executable(simple src/tools/simple.cpp)
    target_link_libraries(simple oink)

    add_executable(rngame src/tools/rngame.cpp)

    add_executable(stgame src/tools/stgame.cpp)
    target_link_libraries(stgame oink)

    add_executable(counter_rr src/tools/counter_rr.cpp)
    target_link_libraries(counter_rr oink)

    add_executable(counter_dp src/tools/counter_dp.cpp)
    target_link_libraries(counter_dp oink)

    add_executable(counter_m src/tools/counter_m.cpp)
    target_link_libraries(counter_m oink)

    add_executable(counter_core src/tools/counter_core.cpp)
    target_link_libraries(counter_core oink)

    add_executable(counter_rob src/tools/counter_rob.cpp)
    target_link_libraries(counter_rob oink)

    add_executable(tc src/tools/tc.cpp)
    target_link_libraries(tc oink)

    add_executable(tc+ src/tools/tc+.cpp)
    target_link_libraries(tc+ oink)

    add_executable(counter_qpt src/tools/counter_qpt.cpp)
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/oink.pc.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/oink.pc" @ONLY)

install(TARGETS oink DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(TARGETS solve DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES src/oink.hpp src/game.hpp src/uintqueue.hpp DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/oink.pc" DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

