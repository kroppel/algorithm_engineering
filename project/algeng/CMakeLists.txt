cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(algeng LANGUAGES CXX)

# search for OpenMP, REQUIRED -> fatal error if not found
find_package(OpenMP REQUIRED)

# add a library to the project using the specified source files
add_library(algeng STATIC include/quicksort.h src/quicksort.cpp)

# PUBLIC -> targets that link to this target get that include directory
# src directory is made PUBLIC because header file needs to include src file
# to access the implementation of a function template
target_include_directories(algeng PUBLIC include PUBLIC src)

# link to the imported target provided by the FindOpenMP module
target_link_libraries(algeng PUBLIC OpenMP::OpenMP_CXX)