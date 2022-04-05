#include "catch.hpp"
#include "quicksort.h"

TEST_CASE("Quicksort: vector of type <int>", "[correctness]") {
    std::vector<int> v, w;
    int values[5] {2, 4, 3, 1, 0};

    for (const int value:values) {
        v.push_back(value);
    }
    for (int i=0; i<5; ++i) {
        w.push_back(i);
    }
    quicksort_parallelized(v, 0, v.size()-1);
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
    quicksort_parallelized(v, 0, v.size()-1);
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
    quicksort_parallelized(v, 0, v.size()-1);
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