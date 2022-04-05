#ifndef ALG_ENG_FUNCTIONS
#define ALG_ENG_FUNCTIONS
#include "quicksort.h"
#endif //ALG_ENG_FUNCTIONS
#include <omp.h>
#include <iostream>
#include <random>       // Header for random number generation
#include <vector>

// Minimum number of vector elements for a vector to be processed by multiple threads
const int MINIMUM_VECTOR_ELEMENT_NUMBER = 100000;

using namespace std;

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

// Retrieve the element from v that has index k in sorted vector v'
// Returns: Element v'[k]
template <typename T>
T quickselect(vector<T>& v, int l_bound, int u_bound, int k) {
    if (l_bound == u_bound) {
        return v.at(l_bound);
    }
    int p = u_bound;
    p = partition(v, l_bound, u_bound, p);
    if (p == k) {
        return v.at(p);
    }
    else if (p < k) {
        return quickselect(v, p+1, u_bound, k);
    }
    else {
        return quickselect(v, l_bound, p-1, k);
    }
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

