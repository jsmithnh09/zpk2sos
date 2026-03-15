#------------------------------------------
# LIBZPK2SOS generation; static for tests.
#------------------------------------------
add_library(zpk2sos SHARED zpk2sos.c)
target_compile_features(zpk2sos PUBLIC cxx_std_11)
if (WIN32 AND MSVC)
	target_compile_options(zpk2sos PUBLIC /W4)
else()
	target_compile_options(zpk2sos PUBLIC -Wall -O2)
	target_link_libraries(zpk2sos m)
endif()
target_include_directories(zpk2sos PUBLIC 
	${CMAKE_CURRENT_SOURCE_DIR}/../include)

# setting properties for Python testing (toplevel already set
# C11 standard requirement.)
set_target_properties(zpk2sos PROPERTIES
	POSITION_INDEPENDENT_CODE ON
	)

