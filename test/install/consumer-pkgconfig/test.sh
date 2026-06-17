#!/bin/sh
set -e

echo "Testing static linking with pkg-config --static..."

echo "CFLAGS:"
pkg-config --cflags --static oink

echo "LIBS:"
pkg-config --libs --static oink

echo "Compiling with static link..."
# Oink's headers require C++17; the .pc deliberately does not force -std.
c++ -std=c++17 main.cpp $(pkg-config --cflags --static oink) $(pkg-config --libs --static oink) -o test

echo "Running binary..."
./test

echo "pkg-config consumer test PASSED"
