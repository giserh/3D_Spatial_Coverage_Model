cmake_minimum_required(VERSION 3.14)
project(coverage)

set(CMAKE_CXX_STANDARD 14)

add_executable(alg1 alg1.cpp )
add_executable(alg3 alg3.cpp)

add_executable(main main.cpp database.cpp database.h)