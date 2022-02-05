cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(my_test LANGUAGES CXX)

message(${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(my_test main.cpp tests.cpp)
# Add compile option to avoid "too many sections" error when including catch.hpp
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-Wa,-mbig-obj COMPILER_SUPPORT)
if(${COMPILER_SUPPORT} EQUAL 1)
    target_compile_options(my_test PUBLIC -Wa,-mbig-obj)
endif()

target_link_libraries(my_test PRIVATE algeng PRIVATE catch)