#!/bin/bash
# Simple TDD script: build and run tests on changes

mkdir -p build
cd build
cmake ..
make
ctest --output-on-failure