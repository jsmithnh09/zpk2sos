enable_testing()
cmake_minimum_required(VERSION 3.21)
project(zpk2sos
	DESCRIPTION "Z-plane ZPG to SOS pairing."
	VERSION 0.0.1
	LANGUAGES C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)

# Python + pytest for testing
find_package(Python3 REQUIRED COMPONENTS Interpreter)

enable_testing()

add_test(
    NAME python_tests
    COMMAND ${Python3_EXECUTABLE} -m pytest ${CMAKE_SOURCE_DIR}/test
)

# Ensure Python can import the built shared library
set_tests_properties(python_tests PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src"
)
