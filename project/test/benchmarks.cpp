#include "catch.hpp"
#include "quicksort.hpp"
#include <random>
#include <iomanip>
#include <algorithm>

using namespace algeng;

TEST_CASE("Partition_fetch_add: vector of type <int> with 200000000 elements", "[performance+correctness]") {
    double start_time, runtime1, runtime2, runtime3;
    const size_t size = 1000000000;
    std::vector<int> v1, v2, v3;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> d(0, size);
    // Fill vectors with random numbers between 0 and size - 1
    std::cout << "Create vectors with " << size << " random elements.\n";
    for (int i = 0; i < size; ++i) {
    v1.push_back(d(gen));
    }
    std::cout << "Vector created.\n";

    v2 = v1;
    v3 = v1;

    int pivot = v1.at(size-1);

    start_time = omp_get_wtime();
    std::partition(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;});
    runtime1 = omp_get_wtime() - start_time;

    start_time = omp_get_wtime();
    partition_fetch_add(v2, 0, v2.size()-1, v2.size()-1, 16*2048, 2);
    runtime2 = omp_get_wtime() - start_time;

    start_time = omp_get_wtime();
    partition_pivot(v3, 0, v3.size()-1, pivot);
    runtime3 = omp_get_wtime() - start_time;

    REQUIRE(std::is_partitioned(v2.begin(), v2.end(), [&pivot](int x) {return x<=pivot;}));
    REQUIRE(std::is_partitioned(v3.begin(), v3.end(), [&pivot](int x) {return x<=pivot;}));


    std::cout << "std::partition with runtime: " << std::setprecision(6) << runtime1 << "s\n";
    std::cout << "partition_fetch_add with runtime: " << std::setprecision(6) << runtime2 << "s\n";
    std::cout << "partition_pivot with runtime: " << std::setprecision(6) << runtime3 << "s\n";
    std::cout << "and MINIMUM_VECTOR_ELEMENT_NUMBER=" << MINIMUM_VECTOR_ELEMENT_NUMBER << "\n";
}

TEST_CASE("Quicksort_parallel: vector of type <int> with 200000000 elements", "[performance+correctness]") {
    double start_time, runtime1, runtime2, runtime3;
    const size_t size = 200000000;
    std::vector<int> v1, v2, v3;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> d(0, size);
    // Fill vectors with random numbers between 0 and size - 1
    std::cout << "Create vectors with " << size << " random elements.\n";
    for (int i = 0; i < size; ++i) {
        v1.push_back(d(gen));
    }
    std::cout << "Vector created.\n";

    v2 = v1;
    v3 = v1;

    int pivot = v1.at(size-1);

    start_time = omp_get_wtime();
    std::partition(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;});
    runtime1 = omp_get_wtime() - start_time;

    start_time = omp_get_wtime();
    partition_fetch_add(v2, 0, v2.size()-1, v2.size()-1, 16*2048, 2);
    runtime2 = omp_get_wtime() - start_time;

    start_time = omp_get_wtime();
    partition_pivot(v3, 0, v3.size()-1, pivot);
    runtime3 = omp_get_wtime() - start_time;

    REQUIRE(std::is_partitioned(v2.begin(), v2.end(), [&pivot](int x) {return x<=pivot;}));
    REQUIRE(std::is_partitioned(v3.begin(), v3.end(), [&pivot](int x) {return x<=pivot;}));


    std::cout << "std::partition with runtime: " << std::setprecision(6) << runtime1 << "s\n";
    std::cout << "partition_fetch_add with runtime: " << std::setprecision(6) << runtime2 << "s\n";
    std::cout << "partition_pivot with runtime: " << std::setprecision(6) << runtime3 << "s\n";
    std::cout << "and MINIMUM_VECTOR_ELEMENT_NUMBER=" << MINIMUM_VECTOR_ELEMENT_NUMBER << "\n";
}