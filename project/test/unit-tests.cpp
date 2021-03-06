#include "catch.hpp"
#include "quicksort.hpp"
#include <random>
#include <iomanip>
#include <algorithm>

using namespace algeng;

TEST_CASE("Partition_fetch_add: vector size smaller < block_size, single-threaded", "[correctness]") {
    std::vector<int> v1;
    int values[9] = {5,6,3,2,7,1,0,0,0};

    for (const int value:values) {
        v1.push_back(value);
    }

    int pivot = v1.back();

    int p = partition_fetch_add(v1, 0, v1.size()-1, 8, 1);

    REQUIRE(std::is_partitioned(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;}));
    REQUIRE(p == 2);
}


TEST_CASE("Partition_fetch_add: vector size between block_size and 2*block_size, single-threaded", "[correctness]") {
    std::vector<int> v1;
    int values[9] = {5,6,3,2,7,1,1,9,4};

    for (const int value:values) {
        v1.push_back(value);
    }

    int pivot = v1.back();

    int p = partition_fetch_add(v1, 0, v1.size()-1, 8, 1);

    REQUIRE(std::is_partitioned(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;}));
    REQUIRE(p == 4);
}

TEST_CASE("Partition_fetch_add: sorted vector, single-threaded", "[correctness]") {
    std::vector<int> v1;
    int values[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

    for (const int value:values) {
        v1.push_back(value);
    }

    int pivot = v1.back();

    int p = partition_fetch_add(v1, 0, v1.size()-1, 128, 1);

    REQUIRE(std::is_partitioned(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;}));
    REQUIRE(p == 15);
}

TEST_CASE("Partition_fetch_add: reversed sorted vector, single-threaded", "[correctness]") {
    std::vector<int> v1;
    int values[16] = {16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};

    for (const int value:values) {
        v1.push_back(value);
    }

    int pivot = v1.back();

    int p = partition_fetch_add(v1, 0, v1.size()-1, 8, 1);

    REQUIRE(std::is_partitioned(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;}));
    REQUIRE(p == 0);
}

TEST_CASE("Partition_fetch_add: vector size smaller < block_size, multithreaded", "[correctness]") {
    std::vector<int> v1;
    int values[9] = {5,6,3,2,7,1,0,0,0};

    for (const int value:values) {
        v1.push_back(value);
    }

    int pivot = v1.back();

    int p = partition_fetch_add(v1, 0, v1.size()-1, 8, 2);

    REQUIRE(std::is_partitioned(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;}));
    REQUIRE(p == 2);
}


TEST_CASE("Partition_fetch_add: vector size between block_size and 2*block_size, multi-threaded", "[correctness]") {
    std::vector<int> v1;
    int values[9] = {5,6,3,2,7,1,1,9,4};

    for (const int value:values) {
        v1.push_back(value);
    }

    int pivot = v1.back();

    int p = partition_fetch_add(v1, 0, v1.size()-1, 8, 2);

    REQUIRE(std::is_partitioned(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;}));
    REQUIRE(p == 4);
}

TEST_CASE("Partition_fetch_add: sorted vector, multi-threaded", "[correctness]") {
    std::vector<int> v1;
    int values[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

    for (const int value:values) {
        v1.push_back(value);
    }

    int pivot = v1.back();

    int p = partition_fetch_add(v1, 0, v1.size()-1, 4, 2);

    REQUIRE(std::is_partitioned(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;}));
    REQUIRE(p == 15);
}

TEST_CASE("Partition_fetch_add: reversed sorted vector, multi-threaded", "[correctness]") {
    std::vector<int> v1;

    int values[16] = {16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};

    for (const int value:values) {
        v1.push_back(value);
    }

    int pivot = v1.back();

    int p = partition_fetch_add(v1, 0, v1.size()-1, 8, 2);

    REQUIRE(std::is_partitioned(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;}));
    REQUIRE(p == 0);
}

TEST_CASE("Partition_fetch_add: equal vector, multi-threaded", "[correctness]") {
    std::vector<int> v1;
    int values[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

    for (const int value:values) {
        v1.push_back(value);
    }

    int pivot = v1.back();

    int p = partition_fetch_add(v1, 0, v1.size()-1, 8, 2);

    REQUIRE(std::is_partitioned(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;}));
    REQUIRE(p == 15);
}

TEST_CASE("Quicksort: vector of type <int>", "[correctness]") {
    std::vector<int> v, w;
    int values[5] {2,4,3,1,0};

    for (const int value:values) {
        v.push_back(value);
    }
    for (int i=0; i<5; ++i) {
        w.push_back(i);
    }
    quicksort_parallel_wrapper(v, 4096, 2);
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
    quicksort_parallel_wrapper(v, 4096, 2);
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
    quicksort_parallel_wrapper(v, 4096, 2);
    REQUIRE(v==w);
}

TEST_CASE("Quickselect: vector of type <int>", "[correctness]") {
    std::vector<int> v, w;
    int k = 3;
    int values[5] {2, 4, 3, 1, 0};


    for (const int value:values) {
        v.push_back(value);
    }

    int vk = quickselect_parallel_wrapper(v, k, 4096, 2);
    REQUIRE(vk==3);
}

TEST_CASE("Quickselect: vector of type <double>", "[correctness]") {
    std::vector<double> v, w;
    int k = 3;
    double values[5] {1.67, 2.38, 0.4, 1.11, 5.9};

    for (const double value:values) {
        v.push_back(value);
    }

    double vk = quickselect_parallel_wrapper(v, k, 4096, 2);
    REQUIRE(vk==2.38);
}

TEST_CASE("Quickselect: vector of type <char>", "[correctness]") {
    std::vector<char> v, w;
    int k = 3;
    char values[5] {'e', 'x', 'a', 'd', 'h'};

    for (const char value:values) {
        v.push_back(value);
    }

    char vk = quickselect_parallel_wrapper(v, k, 4096, 2);
    REQUIRE(vk=='h');
}