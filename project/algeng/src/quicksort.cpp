#ifndef ALG_ENG_FUNCTIONS
#define ALG_ENG_FUNCTIONS
#include "quicksort.h"
#endif //ALG_ENG_FUNCTIONS
#include <omp.h>
#include <iostream>
#include <vector>
#include <atomic>
#include <functional>
#include <iterator>
#include <algorithm>

// Minimum number of vector elements for a vector to be processed by multiple threads
const int MINIMUM_VECTOR_ELEMENT_NUMBER = 100000;

// Partition (sub-)vector v[l_bound:u_bound] on element with index p
// Returns: index of pivot element after partitioning
template <typename T>
int partition(std::vector<T>& v, int l_bound, int u_bound, int p) {
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
int partition_pivot(std::vector<T>& v, int l_bound, int u_bound, T pivot) {
    T buffer;
    int i = l_bound;
    int j = u_bound;

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
    int return_value = (v.at(i)<= pivot) ? i: i-1;
    return return_value;
}

// Partition (sub-)vector v[l_bound:u_bound] on element with index p
// Returns: index of pivot element after partitioning
template<typename T>
int partition_fetch_add(std::vector<T>& v, int size, int p, int block_size) {
    int num_blocks = size / block_size;
    T pivot = v.at(p);
    int number_of_threads = 2;
    int return_value;

    // Atomic indices for vector access
    std::atomic<int> i(0);
    std::atomic<int> j(0);
    std::atomic<int> k(num_blocks-1);

    std::vector<int> clean_up_left(number_of_threads, -1);
    std::vector<int> clean_up_right(number_of_threads, -1);
    std::atomic<int> cul(0);
    std::atomic<int> cur(0);

#pragma omp parallel num_threads(number_of_threads) shared(num_blocks, pivot, i, j, k)
    {
        bool fetch_left = true;
        bool fetch_right = true;

        int l, r;

        std::cout << "Call Function\n";

        while (std::atomic_fetch_add(&i, 1) < num_blocks) {
            if (fetch_left) {
                l = std::atomic_fetch_add(&j, 1);
                fetch_left = false;
            } else if (fetch_right) {
                r = std::atomic_fetch_add(&k, -1);
                fetch_right = false;
            }

            if (!(fetch_left || fetch_right)) {
                // block offset
                int a = 0;
                int b = 0;
                // global indices of first block element
                int v_a = l * block_size;
                int v_b = r * block_size;
                T buffer;
                bool swap_left = false;
                bool swap_right = false;

                while (a < block_size && b < block_size) {
                    while (a < block_size) {
                        if (v.at(v_a + a) > pivot) {
                            swap_left = true;
                            break;
                        }
                        a++;
                    }
                    while (b < block_size) {
                        if (v.at(v_b + b) <= pivot) {
                            swap_right = true;
                            break;
                        }
                        b++;
                    }
                    if (swap_left && swap_right) {
                        buffer = v.at(v_a + a);
                        v.at(v_a + a) = v.at(v_b + b);
                        v.at(v_b + b) = buffer;
                        swap_left = false;
                        swap_right = false;
                    }
                }
                // Check which blocks where processed completely
                if (a == block_size)
                    fetch_left = true;
                if (b == block_size)
                    fetch_right = true;
            }
        }
        /* Clean up preparation:
         *  if fetch_left = true and fetch_right = true -> all blocks were processed completely -> no clean up needed
         *  if fetch_left = false -> left block needs to be processed further
         *  if fetch_right = false -> right block needs ...
        */
        if (!fetch_left) {
            clean_up_left.at(std::atomic_fetch_add(&cul, 1)) = l * block_size;
        }
        if (!fetch_right) {
            clean_up_right.at(std::atomic_fetch_add(&cur, 1)) = r * block_size;
        }
    }
    std::cout << "First Phase finished\n";
    /*
     * Clean up 1st step: (sorting maybe two threads)
     *  Swap elements between remaining left and right blocks
     *  -> everything left of j or everything right of k is now correctly partitioned
     */
    std::sort(clean_up_left.data(), clean_up_left.data()+cul.load());
    std::sort(clean_up_right.data(), clean_up_right.data()+cur.load(), [](int a, int b) {
        return a > b;});

    int a_block = 0;
    int b_block = 0;

    T buffer;
    bool swap_left = false;
    bool swap_right = false;

    while (a_block < cul.load() && b_block < cur.load()) {
        int a = 0;
        int b = 0;

        while (a < block_size && b < block_size) {
            while (a < block_size) {
                if (v.at(clean_up_left.at(a_block) + a) > pivot) {
                    swap_left = true;
                    break;
                }
                a++;
            }
            while (b < block_size) {
                if (v.at(clean_up_left.at(b_block) + b) <= pivot) {
                    swap_right = true;
                    break;
                }
                b++;
            }
            if (swap_left && swap_right) {
                buffer = v.at(clean_up_left.at(a_block) + a);
                v.at(clean_up_left.at(a_block) + a) = v.at(clean_up_left.at(b_block) + b);
                v.at(clean_up_left.at(b_block) + b) = buffer;
                swap_left = false;
                swap_right = false;
            }
        }
        if (a == block_size)
            a_block++;
        if (b == block_size)
            b_block++;
    }
    std::cout << "Cleanup first step finished\n";

    /*
    * Clean up 2nd step:
    * partition remainder between L and R (maybe do in parallel with 1st step)
    */

    int partition_index = partition_pivot(v, j.load() * block_size, k.load() * block_size + block_size - 1, pivot);

    std::cout << "Cleanup second step finished\n";

    if (a_block == cul.load()) {
        int l = partition_index;

        /*
         * Cleanup 3rd step:
         * Swap elements > pivot starting from pivot border l
         * with elements < pivot from the remaining blocks
         */

        swap_left = false;
        swap_right = false;
        int b = 0;

        while (b_block < cur.load() && (clean_up_right.at(b_block)+block_size-1-b) > l) {
            while (b < block_size) {
                if (v.at(clean_up_right.at(b_block)+block_size-1-b) <= pivot) {
                    swap_right = true;
                    break;
                }
                b++;
            }
            while (l < clean_up_right.at(b_block)+block_size-1-b) {
                if (v.at(l) > pivot) {
                    swap_left = true;
                    break;
                }
                l++;
            }
            if (swap_right && swap_left) {
                buffer = v.at(clean_up_right.at(b_block)+block_size-1-b);
                v.at(clean_up_right.at(b_block)+block_size-1-b) = v.at(l);
                v.at(l) = buffer;
                swap_left = false;
                swap_right = false;
            }
            if (b == block_size) {
                b = 0;
                b_block++;
            }
        }
        return_value = (v.at(l)<=pivot) ? l: l-1;
    }
    else {
        /*
         * Clean up 2nd step:
         * partition remainder between L and R (maybe do in parallel with 1st step)
         */

        int r = partition_index-1;

        /*
         * Cleanup 3rd step:
         * Swap elements > pivot starting from pivot border l
         * with elements < pivot from the remaining blocks
         */

        swap_left = false;
        swap_right = false;
        int a = 0;

        while (a_block < cul.load() && (clean_up_left.at(a_block)+a) < r) {
            while (a < block_size) {
                if (v.at(clean_up_left.at(a_block)+a) <= pivot) {
                    swap_right = true;
                    break;
                }
                a++;
            }
            while (r < clean_up_left.at(a_block)+a) {
                if (v.at(r) > pivot) {
                    swap_left = true;
                    break;
                }
                r++;
            }
            if (swap_right && swap_left) {
                buffer = v.at(clean_up_left.at(a_block)+a);
                v.at(clean_up_left.at(a_block)+a) = v.at(r);
                v.at(r) = buffer;
                swap_left = false;
                swap_right = false;
            }
            if (a == block_size) {
                a = 0;
                a_block++;
            }
        }
        return_value = (v.at(r)<=pivot) ? r: r-1;
    }
    std::cout << "Cleanup third step finished\n";
    return return_value;
}

/*// Partition (sub-)vector v[l_bound:u_bound] on element with index p
// Returns: index of pivot element after partitioning
template <class It>
using T = typename std::iterator_traits<It>::value_type;

template<class It, class Compare = std::less<T<It>>>
int partition_strided(It start, It end, Compare cmp = Compare{}) {
    auto const size = std::distance(start, end);
    int buffer[2*omp_get_num_threads()];
    T<It> buffer_left;
    T<It> buffer_right;

    std::atomic<int> i(0);
    std::atomic<int> j(0);
    std::atomic<int> k(size-1);
    std::atomic<int> b_fetch(0);
    std::atomic<int> b_store(0);
    std::atomic<int> phase1_synch(0);
    std::atomic<int> phase2_synch(0);
    std::atomic<int> phase3_synch(0);


    int l, r;
    bool swap_elements = false;

#pragma omp parallel num_threads(2) shared(v, i, j, k, buffer, b_fetch, b_store, phase1_synch, phase2_synch, phase3_synch) private(l, r, buffer_left, buffer_right) firstprivate(pivot, swap_elements, size)
    {
        std::cout << omp_get_num_threads() << "\n";
        while (int t = atomic_fetch_add(&i, 1) < size) {
            //std::cout << "Thread " << omp_get_thread_num() << "\n";
            if (!swap_elements) {
                l = atomic_fetch_add(&j, 1);
                buffer_left = v.at(l);
                //std::cout << "Buffer left: " << (int) buffer_left << "\n";
                if (!cmp(buffer_left)) {
                    swap_elements = true;
                }
            } else {
                r = atomic_fetch_add(&k, -1);
                buffer_right = v.at(r);
                if (cmp(buffer_right)) {
                    // -> each index > k is guaranteed to hold elements > pivot, as each element <= pivot
                    // gets switched with and element lower than j at some point
                    v.at(l) = buffer_right;
                    v.at(r) = buffer_left;
                    swap_elements = false;
                    //std::cout << "Swapping: " << (int) buffer_left << " and " << (int) buffer_right << "\n";
                }
            }
        }
//#pragma omp barrier
        atomic_fetch_add(&phase1_synch, 1);
        while(phase1_synch.load() < omp_get_thread_num()) {

        }
        // before this step j holds the number of left side elements < pivot or that are > pivot but did
        // not find the match to get switched
        if (swap_elements) {
            atomic_fetch_add(&j,-1);
        }
        // after this step j holds the number of left side elements < pivot, as each process with swap_elements=true
        // decrements j (swap_elements=true means that the process found an index j with v[j] < pivot, but no match to switch)

//#pragma omp barrier
        atomic_fetch_add(&phase2_synch, 1);
        while(phase2_synch.load() < omp_get_thread_num()) {

        }


        if (swap_elements) {
            if (l<j.load()) {
                r = atomic_fetch_add(&k,-1);
                if (v.at(r) < pivot && r > j.load()) {
                    buffer[atomic_fetch_add(&b_fetch,1)] = r;
                }
            }
            // processes with l >= j do not need to swap, as there are not enough elements to swap and their left index l
            // is to the right of the cutting point (pivot point)
            else {
                swap_elements = false;
            }
        }

//#pragma omp barrier

        atomic_fetch_add(&phase3_synch, 1);
        while(phase3_synch.load() < omp_get_thread_num()) {

        }

        if (swap_elements) {
            r = buffer[atomic_fetch_add(&b_store,1)];
            buffer_right = v.at(r);
            v.at(l) = buffer_right;
            v.at(r) = v.at(l);
        }
    }

    return i.load();
}
 */

// Retrieve the element from v that has index k in sorted vector v'
// Returns: Element v'[k]
template <typename T>
T quickselect(std::vector<T>& v, int l_bound, int u_bound, int k) {
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
void quicksort(std::vector<T>& v, int l_bound, int u_bound) {
    if (u_bound > l_bound) {
        int p = partition(v, l_bound, u_bound, u_bound);

        quicksort(v, l_bound, p - 1);
        quicksort(v, p + 1, u_bound);
    }
}

template <typename T>
void quicksort_parallel(std::vector<T>& v, int l_bound, int u_bound) {
    if (u_bound > l_bound) {
        int p = partition(v, l_bound, u_bound, u_bound);

#pragma omp parallel sections
        {
#pragma omp section
            {
                //std::cout << "Thread " << omp_get_thread_num() << " starting his task" << "\n";
                if (p - l_bound > MINIMUM_VECTOR_ELEMENT_NUMBER)
                    quicksort_parallel(v, l_bound, p - 1);
                else
                    quicksort(v, l_bound, p - 1);
            }
#pragma omp section
            {
                //std::cout << "Thread " << omp_get_thread_num() << " starting his task" << "\n";
                if (u_bound - p > MINIMUM_VECTOR_ELEMENT_NUMBER)
                    quicksort_parallel(v, p + 1, u_bound);
                else
                    quicksort(v, p + 1, u_bound);
            }
        }
    }
}

