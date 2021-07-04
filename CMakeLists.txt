# Check if we are solo or as part of a larger project
# If we're a subproject, disable testing/benchmark tools
if(DEFINED PROJECT_NAME)
  set(subproject ON)
else()
  set(subproject OFF)
endif()

cmake_minimum_required(VERSION 3.13)

project(oink
    VERSION 1.0.0
    DESCRIPTION "Oink, an implementation of modern parity game solvers"
    HOMEPAGE_URL "https://github.com/trolando/oink"
    LANGUAGES C CXX
)

# Default build type RelWithDebInfo
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Print some info
message(STATUS "CMake build configuration for Oink ${PROJECT_VERSION}")
message(STATUS "OS: ${CMAKE_SYSTEM_NAME}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")

# require at least gcc 4.9
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9)
    message(FATAL_ERROR "GCC version must be at least 4.9! (found ${CMAKE_CXX_COMPILER_VERSION}")
endif()

# Make shared or static library?
option(BUILD_SHARED_LIBS "Build a shared library" OFF)
if(DEFINED oink_SHARED_LIBS)
    set(BUILD_SHARED_LIBS "${oink_SHARED_LIBS}")
endif()

# Required package
find_package(Boost REQUIRED filesystem iostreams regex system random)
include(GNUInstallDirs)

# Main library target
add_library(oink)
add_library(oink::oink ALIAS oink)
set_target_properties(oink PROPERTIES VERSION ${oink_VERSION} SOVERSION ${oink_VERSION_MAJOR})
target_compile_features(oink PUBLIC c_std_11 cxx_std_11)
target_compile_options(oink PRIVATE -Wall -Wextra)
target_compile_options(oink PRIVATE "$<$<CONFIG:Release,RelWithDebInfo>:-march=native>")
target_include_directories(oink
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/src
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src ${Boost_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:include/oink>
)
target_link_libraries(oink pthread)

# set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

target_sources(oink
  PRIVATE
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
    src/dftl.cpp
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src/oink.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src/error.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src/game.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src/bitset.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src/uintqueue.hpp>
)

get_target_property(oink_PUBLIC_HEADERS oink INTERFACE_SOURCES)
set_target_properties(oink PROPERTIES PUBLIC_HEADER "${oink_PUBLIC_HEADERS}")


option(OINK_BUILD_TOOLS "Build basic Oink tools" ON)
option(OINK_BUILD_EXTRA_TOOLS "Build extra Oink tools" OFF)

# Common target properties
macro(set_target_props NAME)
    target_compile_features(${NAME} PUBLIC c_std_11 cxx_std_11)
    target_compile_options(${NAME} PRIVATE -Wall -Wextra)
    target_include_directories(${NAME} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/src ${Boost_INCLUDE_DIRS})
    target_link_libraries(${NAME} oink::oink pthread)
endmacro(set_target_props)


if(OINK_BUILD_TOOLS)
    add_executable(solve)
    target_sources(solve PRIVATE src/tools/solve.cpp src/tools/getrss.c)
    set_target_properties(solve PROPERTIES OUTPUT_NAME oink)
    set_target_props(solve)
    target_link_libraries(solve Boost::iostreams Boost::random)

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
    target_link_libraries(counter_rr oink)

    add_executable(counter_dp src/tools/counter_dp.cpp)
    set_target_props(counter_dp) 
    target_link_libraries(counter_dp oink)

    add_executable(counter_m src/tools/counter_m.cpp)
    set_target_props(counter_m) 
    target_link_libraries(counter_m oink)

    add_executable(counter_core src/tools/counter_core.cpp)
    set_target_props(counter_core) 
    target_link_libraries(counter_core oink)

    add_executable(counter_rob src/tools/counter_rob.cpp)
    set_target_props(counter_rob) 
    target_link_libraries(counter_rob oink)

    add_executable(tc src/tools/tc.cpp)
    set_target_props(tc) 
    target_link_libraries(tc oink)

    add_executable(tc+ src/tools/tc+.cpp)
    set_target_props(tc+) 
    target_link_libraries(tc+ oink)

    add_executable(counter_qpt src/tools/counter_qpt.cpp)
    set_target_props(counter_qpt) 

    add_executable(counter_ortl src/tools/counter_ortl.cpp)
    set_target_props(counter_ortl)
    target_link_libraries(counter_ortl oink)

    add_executable(counter_dftl src/tools/counter_dftl.cpp)
    set_target_props(counter_dftl)
    target_link_libraries(counter_dftl oink)
endif()


if(NOT subproject)
    enable_testing()
    set(CTEST_TESTING_TIMEOUT 10) 

    # test fixpoint algorithms
    add_test(NAME TestSolverFPIpar COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpi -w 0)
    add_test(NAME TestSolverFPIseq COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpi -w -1)
    add_test(NAME TestSolverFPJ COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpj)
    add_test(NAME TestSolverFPJG COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --fpjg)
    # test parallel strategy iteration
    add_test(NAME TestSolverPSIs COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --psi -w -1)
    add_test(NAME TestSolverPSIp COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --psi -w 0)
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
    add_test(NAME TestSolverRRDP COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rrdp)
    # test TL and variations
    add_test(NAME TestSolverTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --tl)
    add_test(NAME TestSolverRTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --rtl)
    add_test(NAME TestSolverORTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ortl)
    add_test(NAME TestSolverPTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --ptl)
    add_test(NAME TestSolverSPPTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --spptl)
    add_test(NAME TestSolverDTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --dtl)
    add_test(NAME TestSolverIDTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --idtl)
    add_test(NAME TestSolverDFTL COMMAND test_solvers ${CMAKE_CURRENT_SOURCE_DIR}/tests --dftl)
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


if(NOT subproject)
  install(TARGETS oink
    EXPORT OinkConfig
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}"
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}"
  )
  install(TARGETS solve
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  )
  install(EXPORT OinkConfig FILE "OinkConfig.cmake" NAMESPACE oink:: DESTINATION cmake)

  # Setup basic CPack support
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
  set(CMAKE_PROJECT_HOMEPAGE_URL ${PROJECT_URL})
  set(CPACK_PACKAGE_VENDOR "Tom van Dijk")
  set(CPACK_PACKAGE_CONTACT "Tom van Dijk <tom@tvandijk.nl>")
  set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.md")
  set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
  set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.md")
  set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
  set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_MAJOR_VERSION})
  set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_MINOR_VERSION})
  set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_PATCH_VERSION})
  include(CPack)
endif()
