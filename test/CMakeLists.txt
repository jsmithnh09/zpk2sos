# Python + pytest for testing
find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(
    NAME python_tests
    COMMAND ${Python3_EXECUTABLE} -m pytest ${CMAKE_SOURCE_DIR}/test
)

# Ensure Python can import the built shared library
set_tests_properties(python_tests PROPERTIES
    ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src"
)