#include "catch.hpp"
#include "quicksort.h"
#include <random>
#include <iomanip>
#include <algorithm>

TEST_CASE("Partition_fetch_add: vector size smaller < block_size, single-threaded", "[correctness]") {
    std::vector<int> v1, v2;
    int values[9] = {5,6,3,2,7,1};

    for (const int value:values) {
        v1.push_back(value);
    }

    v2 = v1;
    int pivot = v1.back();

    partition(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;});
    partition_fetch_add(v2, v2.size(), v2.size()-1, 8, 1);

    REQUIRE(v1 == v2);
    REQUIRE(std::is_partitioned(v2.begin(), v2.end(), [&pivot](int x) {return x<=pivot;}));
}


TEST_CASE("Partition_fetch_add: vector size between block_size and 2*block_size, single-threaded", "[correctness]") {
    std::vector<int> v1, v2;
    int values[9] = {5,6,3,2,7,1,1,9,4};

    for (const int value:values) {
        v1.push_back(value);
    }

    v2 = v1;
    int pivot = v1.back();

    partition(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;});
    partition_fetch_add(v2, v2.size(), v2.size()-1, 8, 1);

    REQUIRE(v1 == v2);
    REQUIRE(std::is_partitioned(v2.begin(), v2.end(), [&pivot](int x) {return x<=pivot;}));
}

TEST_CASE("Partition_fetch_add: sorted vector, single-threaded", "[correctness]") {
    std::vector<int> v1, v2;
    int values[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

    for (const int value:values) {
        v1.push_back(value);
    }

    v2 = v1;
    int pivot = v1.back();

    partition(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;});
    partition_fetch_add(v2, v2.size(), v2.size()-1, 128, 1);

    REQUIRE(v1 == v2);
    REQUIRE(std::is_partitioned(v2.begin(), v2.end(), [&pivot](int x) {return x<=pivot;}));
}

TEST_CASE("Partition_fetch_add: reversed sorted vector, single-threaded", "[correctness]") {
    std::vector<int> v1, v2;
    int values[16] = {16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};

    for (const int value:values) {
        v1.push_back(value);
    }

    v2 = v1;
    int pivot = v1.back();

    partition(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;});
    partition_fetch_add(v2, v2.size(), v2.size()-1, 8, 1);

    REQUIRE(v1 == v2);
    REQUIRE(std::is_partitioned(v2.begin(), v2.end(), [&pivot](int x) {return x<=pivot;}));
}

TEST_CASE("Partition_fetch_add: vector of type <int> with 200000000 elements", "[correctness]") {
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
    partition(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;});
    runtime1 = omp_get_wtime() - start_time;

    start_time = omp_get_wtime();
    partition_fetch_add(v2, v2.size(), v2.size()-1, 16*2048, 2);
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

TEST_CASE("Quicksort: vector of type <int>", "[correctness]") {
    std::vector<int> v, w;
    int values[5] {2, 4, 3, 1, 0};

    for (const int value:values) {
        v.push_back(value);
    }
    for (int i=0; i<5; ++i) {
        w.push_back(i);
    }
    quicksort_parallel(v, 0, v.size()-1);
    REQUIRE(v==w);
}

TEST_CASE("Quicksort: vector of type <double>", "[correctness]") {
    std::vector<double> v, w;
    double values[5] {1.67, 2.38, 0.4, 1.11, 5.9};
    double sorted_values[5] {0.4, 1.11, 1.67, 2.38, 5.9};

    for (const double value:values) {
        v.push_back(value);
    }
    for (const double value:sorted_values) {
        w.push_back(value);
    }
    quicksort_parallel(v, 0, v.size()-1);
    REQUIRE(v==w);
}

TEST_CASE("Quicksort: vector of type <char>", "[correctness]") {
    std::vector<char> v, w;
    char values[5] {'e', 'x', 'a', 'd', 'h'};
    char sorted_values[5] {'a', 'd', 'e', 'h', 'x'};

    for (const char value:values) {
        v.push_back(value);
    }
    for (const char value:sorted_values) {
        w.push_back(value);
    }
    quicksort_parallel(v, 0, v.size()-1);
    REQUIRE(v==w);
}

TEST_CASE("Quickselect: vector of type <int>", "[correctness]") {
    std::vector<int> v, w;
    int k = 3;
    int values[5] {2, 4, 3, 1, 0};


    for (const int value:values) {
        v.push_back(value);
    }

    int vk = quickselect(v, 0, v.size()-1, k);
    REQUIRE(vk==3);
}

TEST_CASE("Quickselect: vector of type <double>", "[correctness]") {
    std::vector<double> v, w;
    int k = 3;
    double values[5] {1.67, 2.38, 0.4, 1.11, 5.9};

    for (const double value:values) {
        v.push_back(value);
    }

    double vk = quickselect(v, 0, v.size()-1, k);
    REQUIRE(vk==2.38);
}

TEST_CASE("Quickselect: vector of type <char>", "[correctness]") {
    std::vector<char> v, w;
    int k = 3;
    char values[5] {'e', 'x', 'a', 'd', 'h'};

    for (const char value:values) {
        v.push_back(value);
    }

    char vk = quickselect(v, 0, v.size()-1, k);
    REQUIRE(vk=='h');
}