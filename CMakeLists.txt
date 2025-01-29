cmake_minimum_required(VERSION 3.16)
project(
    oxsx
    VERSION 1.4.0
    DESCRIPTION "The 'Oxford' Signal Extraction framework for the SNO+ experiment"
    LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(OXSX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
file(MAKE_DIRECTORY ${OXSX_INCLUDE_DIR})

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(test)