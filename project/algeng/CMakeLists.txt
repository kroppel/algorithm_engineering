cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(algeng LANGUAGES CXX)

# search for OpenMP, REQUIRED -> fatal error if not found
find_package(OpenMP REQUIRED)

# add an interface library -> does not generate build output by itself
add_library(algeng INTERFACE)

# INTERFACE -> include directories will be inherited by targets that link with this library
target_include_directories(algeng INTERFACE include)

# INTERFACE -> each target that links with this library will also link with OpenMP
target_link_libraries(algeng INTERFACE OpenMP::OpenMP_CXX)