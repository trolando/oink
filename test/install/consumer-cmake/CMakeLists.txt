cmake_minimum_required(VERSION 3.14)
project(oink-consumer-test LANGUAGES C CXX)

find_package(oink REQUIRED)

add_executable(test-oink main.cpp)
target_link_libraries(test-oink PRIVATE oink::oink)
