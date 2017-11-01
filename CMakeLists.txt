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

include_directories(src)
include(GNUInstallDirs)

option(BUILD_SHARED_LIBS "Enable/disable creation of shared libraries" OFF)
option(BUILD_STATIC_LIBS "Enable/disable creation of static libraries" ON)
option(BUILD_EXTRA_TOOLS "Enable/disable building of extra tools" OFF)

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
    src/pp.cpp
    src/ppp.cpp
    src/rr.cpp
    src/dp.cpp
    src/rrdp.cpp
    src/apt.cpp
    src/psi.cpp
    src/spm.cpp
    src/mspm.cpp
    src/tspm.cpp
    src/qpt.cpp
)
target_link_libraries(oink pthread)

add_executable(solve src/tools/solve.cpp)
set_target_properties(solve PROPERTIES OUTPUT_NAME oink)
target_link_libraries(solve oink boost_iostreams)

if(BUILD_EXTRA_TOOLS)
    add_executable(nudge src/tools/nudge.cpp)
    target_link_libraries(nudge oink)

    add_executable(dotty src/tools/dotty.cpp)
    target_link_libraries(dotty oink)

    add_executable(verify src/tools/verify.cpp)
    target_link_libraries(verify oink)

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

    add_executable(counter_qpt src/tools/counter_qpt.cpp)
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/oink.pc.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/oink.pc" @ONLY)

install(TARGETS oink DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(TARGETS solve DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES src/oink.hpp src/game.hpp DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/oink.pc" DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

