cmake_minimum_required(VERSION 3.0.0)

project(BLAS)

include_directories(/home/mjdyx/Documents/C++/include)
add_compile_options(-std=c++17)
add_compile_options(-march=native)

add_subdirectory(MatMultMT)
add_subdirectory(TestSet)