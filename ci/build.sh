#!/usr/bin/bash
#
# Intended for running in CI

set -Eeuo pipefail

# groovy one liner.
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

cd "$SCRIPT_DIR/.."
GREEN_COLOR='\033[0;32m'
RESET_COLOR='\033[0m'
CHECKMARK='\u2714'

if ! command -v cmake &> /dev/null
then
	echo "CMake could not be found."
	exit 1
fi
echo "-- CMake installed."
if ! command -v uv &> /dev/null
then
	echo "uv python manager not found."
	exit 1
fi
echo "-- uv installed."

# generate the virtual environment for testing.
uv venv --clear
uv pip install -e .

# Build the files.
uv run cmake --preset test
uv run cmake --build --preset test

# Running CTest.
uv run ctest --preset test -V --output-on-failure

# hurray!
printf "${GREEN_COLOR}${CHECKMARK} Success!${RESET_COLOR}\n"
exit 0
