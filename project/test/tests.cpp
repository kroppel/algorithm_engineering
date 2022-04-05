#include "catch.hpp"
#include "quicksort.h"

TEST_CASE("Vector of type <int>", "[correctness]") {
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

TEST_CASE("Vector of type <double>", "[correctness]") {
    std::vector<int> v, w;
    double values[5] {1.67, 2.38, 0.4, 1.11, 5.9};
    double sorted_values[5] {0.4, 1.11, 1.67, 2.38, 5.9};

    for (const int value:values) {
        v.push_back(value);
    }
    for (const int value:sorted_values) {
        w.push_back(value);
    }
    quicksort_parallelized(v, 0, v.size()-1);
    REQUIRE(v==w);
}

TEST_CASE("Vector of type <Char>", "[correctness]") {
    std::vector<int> v, w;
    char values[5] {'e', 'x', 'a', 'd', 'h'};
    char sorted_values[5] {'a', 'd', 'e', 'h', 'x'};

    for (const int value:values) {
        v.push_back(value);
    }
    for (const int value:sorted_values) {
        w.push_back(value);
    }
    quicksort_parallelized(v, 0, v.size()-1);
    REQUIRE(v==w);
}