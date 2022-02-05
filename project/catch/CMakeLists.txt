cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(catch LANGUAGES CXX)

add_library(catch INTERFACE)
target_include_directories(catch INTERFACE include)