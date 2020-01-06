cmake_minimum_required(VERSION 3.15)
project(Alignment)

set(CMAKE_CXX_STANDARD 98)

add_executable(Alignment main.cpp Pair_HMM.h Pair_HMM.cpp helper_functions.cpp helper_functions.h)