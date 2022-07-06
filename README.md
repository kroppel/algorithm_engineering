# Algorithm Engineering
This repository contains all code written for the course "algorithm engineering".

## Introduction

The project consists of a Header-Only-Library providing functions that perform parallel in-place Partitioning, Quickselect and Quicksort. In case of Partitioning and Quickselect, benchmarks show that a major performance gain can be obtained in comparison to serial sorting for sufficiently high workloads, given that enough compute resources are available. All three implementations were benchmarked against *std::sort* and *__gnu_parallel::sort*.

## Usage

To be able to use the functions inside your C++ code, just include the header that can be found under *project/algeng/include/quicksort.hpp*.

### Partitioning

### Quickselect

### Quicksort

## Further Ressources

For further information on the implementation, look into the specific source code and consult the paper that can be found in the directory *paper/*
