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