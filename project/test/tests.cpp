#include "catch.hpp"
#include "quicksort.h"
#include <random>
#include <iomanip>
#include <algorithm>

using namespace std;

TEST_CASE("Partition_fetch_add: vector of type <int> with 6 elements", "[correctness]") {
    vector<int> v1, v2;
    int values[6] = {5,6,3,2,7,1};

    for (const int value:values) {
        v1.push_back(value);
    }

    v2 = v1;
    int pivot = v1.back();
    cout << "Pivot: " << pivot << "\n";

    partition(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;});

    partition_fetch_add(v2, v2.size(), pivot);

    REQUIRE(v1 == v2);
}

TEST_CASE("Partition_fetch_add: vector of type <int> with 1000000 elements", "[correctness]") {
    double start_time, runtime1, runtime2;
    int p1, p2;
    const size_t size = 10;
    vector<int> v1, v2;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> d(0, size);
    // Fill vectors with random numbers between 0 and size - 1
    cout << "Create vectors with " << size << " random elements.\n";
    for (int i = 0; i < size; ++i) {
        v1.push_back(d(gen));
    }

    v2 = v1;
    int pivot = v1.at(size-1);
    cout << "Pivot: " << pivot << "\n";

    start_time = omp_get_wtime();
    partition(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;});
    runtime1 = omp_get_wtime() - start_time;

    start_time = omp_get_wtime();
    partition_fetch_add(v2, v2.size(), pivot);
    runtime2 = omp_get_wtime() - start_time;

    REQUIRE(v1 == v2);

    cout << "std::partition with runtime: " << setprecision(6) << runtime1 << "s\n";
    cout << "partition_fetch_add with runtime: " << setprecision(6) << runtime2 << "s\n";
    cout << "and MINIMUM_VECTOR_ELEMENT_NUMBER=" << MINIMUM_VECTOR_ELEMENT_NUMBER << "\n";
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