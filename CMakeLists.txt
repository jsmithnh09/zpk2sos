cmake_minimum_required(VERSION 3.21)
project(zpk2sos
	DESCRIPTION "Z-plane ZPG to SOS pairing."
	VERSION 0.0.1
	LANGUAGES C)

enable_testing()
add_subdirectory(src)
add_subdirectory(test)

