cmake_minimum_required(VERSION 3.17)
project(TRDP)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")

add_executable(TRDP trdp.cpp)