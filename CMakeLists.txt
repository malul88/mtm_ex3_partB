cmake_minimum_required(VERSION 3.16)
project(mtm_ex3_partB)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic-errors -Werror")

add_executable(mtm_ex3_partB test_partB.cpp Matrix.cpp Auxiliaries.cpp)