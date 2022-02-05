#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_FAST_COMPILE

#include "catch.hpp"

/*#include <omp.h>
#include <vector>
#include <iomanip>
#include <algorithm>    // Header for std::is_sorted function
#include <random>       // Header for random number generation
#include <iostream>
#include "quicksort.h"

using namespace std;

int main() {
    double start_time, runtime1, runtime2;
    const size_t size = 1000000;
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

    start_time = omp_get_wtime();
    quicksort(v1,0, v1.size() - 1);
    runtime1 = omp_get_wtime() - start_time;

    start_time = omp_get_wtime();
    quicksort_parallelized(v2, 0, v2.size() - 1);
    runtime2 = omp_get_wtime() - start_time;

    // Check if vectors are sorted
    cout << "Vector 1 sorted correctly: " << is_sorted(begin(v1), end(v1)) << "\n";
    cout << "With runtime: " << setprecision(6) << runtime1 << "s\n";
    cout << "Vector 2 sorted correctly: " << is_sorted(begin(v2), end(v2)) << "\n";
    cout << "With runtime: " << setprecision(6) << runtime2 << "s\n";
    cout << "and MINIMUM_VECTOR_ELEMENT_NUMBER=" << MINIMUM_VECTOR_ELEMENT_NUMBER << "\n";
}*/