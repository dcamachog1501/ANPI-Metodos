cmake_minimum_required(VERSION 3.19)
project(Metodos)
include(FindPkgConfig)
pkg_check_modules(GiNac REQUIRED ginac)

set(CMAKE_CXX_STANDARD 14)

add_executable(Metodos main.cpp)

target_link_libraries(Metodos ginac)