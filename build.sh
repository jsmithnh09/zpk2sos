#!/usr/bin/bash

cmake -S . -B cmake-build
cd cmake-build
cmake --build .
ctest
exit 0
