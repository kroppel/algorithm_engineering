# Algorithm Engineering

This repository contains all code written for the course "algorithm engineering".

## About this project

The project consists of a Header-Only-Library providing functions that perform parallel in-place Partitioning, Quickselect and Quicksort. In case of Partitioning and Quickselect, benchmarks show that a major performance gain can be obtained in comparison to serial sorting for sufficiently high workloads, given that enough compute resources are available. All three implementations were benchmarked against *std::sort* and *__gnu_parallel::sort*.

## Usage

To be able to use the functions inside your C++ code, just include the header that can be found under *project/algeng/include/quicksort.hpp*.

- __Partitioning__: The function *partition_fetch_add* partitions a given vector using the last element as the pivot. The user can specify the desired block size and number of threads to be used.
  
- __Quickselect__: The function *quickselect_parallel* selects the element out of the given vector that would be located at index k if the vector was sorted. The function quickselect_parallel_wrapper can be used as an alternative to *quickselect_parallel* if one wants to process the whole vector.

- __Quicksort__: The function *quicksort_parallel* sorts a given vector utilizing *partition_fetch_add* as partitioning algorithm for large vectors. The function quicksort_parallel_wrapper can be used as an alternative to *quicksort_parallel* if one wants to process the whole vector.

## Further Ressources

For further information on the implementation, look into the specific source code and consult the paper that can be found in the directory *paper/*.
For the 5 functions listed above there can be more details found in the Docstrings.
