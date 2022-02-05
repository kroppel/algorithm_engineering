cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(test LANGUAGES CXX)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
            "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
            FORCE)
endif()

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add subdirectories for project and testing
add_subdirectory(algeng)
add_subdirectory(catch)

add_executable(test main.cpp tests.cpp)
# Add compile option to avoid "too many sections" error when including catch.hpp
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-Wa,-mbig-obj COMPILER_SUPPORT)
if(${COMPILER_SUPPORT} EQUAL 1)
    target_compile_options(test PUBLIC -Wa,-mbig-obj)
endif()

target_link_libraries(test PRIVATE algeng PRIVATE catch)