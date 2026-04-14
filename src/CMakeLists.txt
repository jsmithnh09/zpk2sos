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
	$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
)

# setting properties for Python testing (toplevel already set
# C11 standard requirement.)
set_target_properties(zpk2sos PROPERTIES
	POSITION_INDEPENDENT_CODE ON
)
set(ZSOS_LIB_PATH "${CMAKE_BINARY_DIR}/libzpk2sos${CMAKE_SHARED_LIBRARY_PREFIX}")
set_property(GLOBAL PROPERTY ENVIRONMENT "ZSOS_LIB_PATH=${ZSOS_LIB_PATH}")


# adding debug capability here (kinda manual, but that's fine for now.)
# option(ENABLE_DEBUG "Enabling Debug printout for pair tracking" ON)
# if (ENABLE_DEBUG)
#	target_compile_definitions(zpk2sos PRIVATE ZPK_DEBUG)
# endif()
