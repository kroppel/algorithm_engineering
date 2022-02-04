#include <omp.h>
#include <iostream>
#include <random>       // Header for random number generation
#include <vector>
#include <algorithm>    // Header for std::is_sorted function
#include <iomanip>

using namespace std;
// Minimum number of vector elements for a vector to be processed by multiple threads
const int MINIMUM_VECTOR_ELEMENT_NUMBER = 100000;

// Partition (sub-)vector v[l_bound:u_bound] on element with index p
// Returns: index of pivot element after partitioning
template <typename T>
int partition(vector<T>& v, int l_bound, int u_bound, int p) {
    T buffer;
    int i = l_bound;
    int j = u_bound - 1;

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
    }

    return i;
}

template <typename T>
void quicksort(vector<T>& v, int l_bound, int u_bound) {
    if (u_bound > l_bound) {
        int i = partition(v, l_bound, u_bound, u_bound);

        quicksort(v, l_bound, i - 1);
        quicksort(v, i + 1, u_bound);
    }
}

template <typename T>
void quicksort_parallelized(vector<T>& v, int l_bound, int u_bound) {
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
#pragma omp parallel sections
        {
            {
                //cout << "Thread " << omp_get_thread_num() << " starting his task" << "\n";
                if (i - l_bound > MINIMUM_VECTOR_ELEMENT_NUMBER)
                    quicksort_parallelized(v, l_bound, i - 1);
                else
                    quicksort(v, l_bound, i - 1);
            }
#pragma omp section
            {
                //cout << "Thread " << omp_get_thread_num() << " starting his task" << "\n";
                if (u_bound - i > MINIMUM_VECTOR_ELEMENT_NUMBER)
                    quicksort_parallelized(v, i + 1, u_bound);
                else
                    quicksort(v, i + 1, u_bound);
            }
        }
    }
}

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
    quicksort(v1, 0, v1.size() - 1);
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
}

