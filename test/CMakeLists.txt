find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(VENV_DIR "${CMAKE_CURRENT_BINARY_DIR}/venv")
set(VENV_PYTHON "${VENV_DIR}/bin/python")

if(NOT EXISTS "${VENV_DIR}")
	message(STATUS "### Creating Python virtual environment")
	execute_process(
		COMMAND ${Python3_EXECUTABLE} -m venv ${VENV_DIR}
		RESULT_VARIABLE venv_result
		)
	if(NOT venv_result EQUAL 0)
		message(FATAL_ERROR "Failed to create environment.")
	endif()
	message(STATUS "### Installing SciPy/PyTest")
	set(PY_PKG "scipy" "pytest" "numpy")
	foreach(pkg ${PY_PKG})
		# checking if the package is already present in the virtual
		# environment; save some build-time and effort.
		execute_process(COMMAND ${VENV_PYTHON} -m pip show ${pkg}
			RESULT_VARIABLE exist_result
			OUTPUT_QUIET ERROR_QUIET
			)
		if(NOT ${exist_result} EQUAL 0)
			execute_process(COMMAND ${VENV_PYTHON} -m pip install ${pkg}
				RESULT_VARIABLE pip_result
				)
			if(NOT ${pip_result} EQUAL 0)
				message(FATAL_ERROR "Failed to install ${pkg}.")
			endif()
		else()
			message(STATUS, "### ${pkg} installed.")
		endif()
	endforeach()
endif()

# adding the test for Python.
add_test(NAME test_coeffs
	COMMAND ${VENV_PYTHON} -u -m pytest --capture=no "${CMAKE_CURRENT_SOURCE_DIR}")

# using environment variable generator instead of copying directly into tree.
set_tests_properties(test_coeffs PROPERTIES ENVIRONMENT
	"ZSOS_LIB_PATH=$<TARGET_FILE:zpk2sos>")
