cmake_minimum_required(VERSION 3.28)
project(zpk2sos
	DESCRIPTION "Z-plane ZPG to SOS pairing."
	VERSION 0.0.1
	LANGUAGES C)


# enable_testing()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_subdirectory(src)
add_subdirectory(test)
