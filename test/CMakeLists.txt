add_executable(test_ellip test_ellip.c)

target_link_libraries(test_ellip PRIVATE zpk2sos)
target_include_directories(test_ellip PUBLIC 
	${CMAKE_CURRENT_SOURCE_DIR}/../include)

# Setup the test.
add_test(NAME TestEllipSos COMMAND test_ellip)
set_tests_properties(TestEllipSos PROPERTIES
	PASS_REGULAR_EXPRESSION "PASS"
	FAIL_REGULUAR_EXPRESSION "FAIL"
	)

