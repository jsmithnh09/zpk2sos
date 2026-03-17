#!/usr/bin/bash
#
# Intended for running (hopefully) in CI??

set -Eeuo pipefail

# groovy one liner.
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

# some constants
SOURCE_DIR=$SCRIPT_DIR/../
BUILD_DIR=$SCRIPT_DIR/../cmake-build
GREEN_COLOR='\033[0;32m'
RESET_COLOR='\033[0m'
CHECKMARK_ICON='\u2714'

# clean build every time
rm -rf "$BUILD_DIR"

# check if Cmake installed
if ! command -v cmake &> /dev/null
then
	echo "CMake could not be found."
	exit 1
fi
echo "-- CMake installed."

# make the build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# run generation/configuration.
echo "-- Generating..."
cmake -S "$SOURCE_DIR" -B .

# run the build.
echo "-- Building..."
cmake --build .

# Run CTest inside the build.
echo "-- Testing..."
ctest -V --output-on-failure

# hurray!
printf "${GREEN_COLOR}${CHECKMARK_ICON} Success!${RESET_COLOR}\n"
exit 0
