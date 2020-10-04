cmake_minimum_required(VERSION 3.0.0)

project(BLAS)

include_directories(/home/mjdyx/Documents/C++/include)
add_compile_options(-std=c++17)
add_compile_options(-march=native)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../bin/Debug)
else()
	set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../bin/Release)
endif()

add_subdirectory(MatMultMT)
add_subdirectory(TestSet)