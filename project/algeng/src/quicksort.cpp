#ifndef ALG_ENG_FUNCTIONS
#define ALG_ENG_FUNCTIONS
#include "quicksort.h"
#endif //ALG_ENG_FUNCTIONS
#include <omp.h>
#include <iostream>
#include <random>       // Header for random number generation
#include <vector>
#include <atomic>

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

// Partition (sub-)vector v[l_bound:u_bound] on pivot p
// Returns: index of pivot element after partitioning
template <typename T>
int partition_pivot(vector<T>& v, int l_bound, int u_bound, T pivot) {
    T buffer;
    int i = l_bound;
    int j = u_bound - 1;

    if (u_bound > l_bound) {
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
    }

    return i;
}

// Partition (sub-)vector v[l_bound:u_bound] on element with index p
// Returns: index of pivot element after partitioning
template <typename T>
int partition_fetch_add(vector<T>& v, int size, int p) {
    vector<T> buffer(2*omp_get_num_threads());
    atomic<int> buffer_index(0);
    T buffer_left, buffer_right;
    T pivot = v.at(p);

    atomic<int> i(0);
    atomic<int> j(0);
    atomic<int> k(size-1);

    int l, r;
    bool swap = false;

    while (atomic_fetch_add(&i,1) < size) {
        T current_element = v.at(i);
        if (!swap) {
            l = atomic_fetch_add(&j, 1);
            buffer_left = v.at(l);
            if (buffer_left > pivot)
                swap = true;
        }
        else {
            r = atomic_fetch_add(&k, -1);
            buffer_right = v.at(r);
            if (buffer_right <= pivot) {
                v.at(l) = buffer_right;
                v.at(r) = buffer_left;
                swap = false;
            }
        }
    }
    /*if (swap) {
        atomic_fetch_add(&j,-1);
    }
    if (swap) {
        if (l<j) {
            r = atomic_fetch_add(&k,-1);
            if (v.at(r) < pivot && r > j) {
                buffer.at(atomic_fetch_add(&buffer_index,1)) = r;
            }
        }
    }*/

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
        int p = partition(v, l_bound, u_bound, u_bound);

        quicksort(v, l_bound, p - 1);
        quicksort(v, p + 1, u_bound);
    }
}

template <typename T>
void quicksort_parallel(vector<T>& v, int l_bound, int u_bound) {
    if (u_bound > l_bound) {
        int p = partition(v, l_bound, u_bound, u_bound);

#pragma omp parallel sections
        {
#pragma omp section
            {
                //cout << "Thread " << omp_get_thread_num() << " starting his task" << "\n";
                if (p - l_bound > MINIMUM_VECTOR_ELEMENT_NUMBER)
                    quicksort_parallel(v, l_bound, p - 1);
                else
                    quicksort(v, l_bound, p - 1);
            }
#pragma omp section
            {
                //cout << "Thread " << omp_get_thread_num() << " starting his task" << "\n";
                if (u_bound - p > MINIMUM_VECTOR_ELEMENT_NUMBER)
                    quicksort_parallel(v, p + 1, u_bound);
                else
                    quicksort(v, p + 1, u_bound);
            }
        }
    }
}

