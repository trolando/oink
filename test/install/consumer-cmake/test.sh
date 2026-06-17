#!/bin/sh
set -e

echo "Testing CMake find_package(oink)..."

cmake -S . -B build
cmake --build build

echo "Running binary..."
./build/test-oink

echo "CMake consumer test PASSED"
