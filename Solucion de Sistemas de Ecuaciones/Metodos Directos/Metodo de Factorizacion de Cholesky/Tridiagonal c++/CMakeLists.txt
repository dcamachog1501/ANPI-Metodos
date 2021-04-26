cmake_minimum_required(VERSION 3.16)
project(untitled)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 14)

add_executable(untitled main.cpp)
target_link_libraries(untitled ${ARMADILLO_LIBRARIES})