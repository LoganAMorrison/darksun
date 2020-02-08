cmake_minimum_required(VERSION 3.1)
project(DARKSUN)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3")

add_subdirectory(darksun)
add_subdirectory(apps)