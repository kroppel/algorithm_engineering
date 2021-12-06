#include <omp.h>
#include <iostream>
#include <random>       // Header for random number generation
#include <vector>
#include <algorithm>    // Header for std::is_sorted function
#include <iomanip>

using namespace std;
template <typename T>

void quicksort_omp(vector<T>& v, int l_bound, int u_bound) {
    T buffer;
    int i = l_bound;
    int j = u_bound - 1;
    int p = u_bound;

    if (u_bound > l_bound) {
        T pivot = v.at(p);

        while (i < j) {
            // Make sure '>' and '<' operators are defined for given T!
            while ((v.at(i) <= pivot) && (i < u_bound)) {
                i++;
            }
            while ((v.at(j) > pivot) && (j > l_bound)) {
                j--;
            }
            if (i < j) {
                // Switch elements at positions i and j
                buffer = v.at(i);
                v.at(i) = v.at(j);
                v.at(j) = buffer;
            }
        }
        if (v.at(i) > pivot) {
            v.at(p) = v.at(i);
            v.at(i) = pivot;
        }

        quicksort_omp(v, l_bound, i-1);
        quicksort_omp(v, i+1, u_bound);
    }
}

int main() {
    const size_t size = 2000000;
    vector<int> v;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> d(0, size);

    // Fill array with random numbers between 0 and size
    cout << "Create vector with " << size << " random elements.\n";
    for (int i = 0; i < size; ++i) {
        v.push_back(d(gen));
    }

    double start_time = omp_get_wtime();
    quicksort_omp(v, 0, v.size() - 1);
    double runtime = omp_get_wtime() - start_time;

    // Check if array is sorted
    cout << "Vector sorted correctly: " << is_sorted(begin(v), end(v)) << "\n";
    cout << "With runtime: " << setprecision(6) << runtime << "s\n";
}

