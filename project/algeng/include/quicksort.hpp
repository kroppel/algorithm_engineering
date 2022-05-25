#ifndef ALGENG
#define ALGENG
#include <omp.h>
#include <iostream>
#include <vector>
#include <atomic>
#include <functional>
#include <iterator>
#include <algorithm>
#endif //ALGENG

namespace algeng {
    const int MINIMUM_VECTOR_ELEMENT_NUMBER = 100000;
    const int MINIMUM_BLOCK_NUMBER = 250;

    template<typename T>
    inline void move_pivot_to_last_index(std::vector<T> &v, int l_bound, int u_bound) {
        // Compute median of 5 elements that are distributed evenly over vector and swap
        // the element with the element on index 'u_bound'.
        // Worst Case O(N^2) Search for Median, but it is only 5 elements.
        int p1[5] = {l_bound, (l_bound+((u_bound+l_bound)/2))/2, (u_bound+l_bound)/2, (((u_bound+l_bound)/2)+u_bound)/2, u_bound};
        for (int a:p1) {
            int sorted_index=0;
            for (int b:p1) {
                if (v.at(a)>v.at(b))
                    sorted_index++;
            }
            if (sorted_index==2) {
                T buff = v.at(u_bound);
                v.at(u_bound) = v.at(a);
                v.at(a) = buff;
                return;
            }
        }
        T buff = v.at(u_bound);
        v.at(u_bound) = v.at((u_bound+l_bound)/2);
        v.at((u_bound+l_bound)/2) = buff;
        return;
    }

// Partition (sub-)vector v[l_bound:u_bound] on element with index p
// Returns: index of pivot element after partitioning
    template<typename T>
    inline int partition(std::vector<T> &v, int l_bound, int u_bound, int p) {
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
    template<typename T>
    inline int partition_pivot(std::vector<T> &v, int l_bound, int u_bound, T pivot) {
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
        int return_value = (v.at(i) <= pivot) ? i : i - 1;
        return return_value;
    }

    template<typename T>
    inline void insertion_sort(std::vector<T> &v, int l_bound, int u_bound) {
        int j;
        T current;
        for (int i = l_bound; i <= u_bound; i++) { // i=3
            current = v.at(i); //1
            j = i - 1; //j=-1

            while (j >= l_bound && v.at(j) > current) {
                v.at(j+1) = v.at(j);
                j = j - 1;//j=-1
            }

            v.at(j+1) = current;
        }
    }

    template<typename T>
    inline void quicksort(std::vector<T> &v, int l_bound, int u_bound) {
        if (u_bound > l_bound) {
            if (u_bound - l_bound + 1 < 64) { // use insertion sort for small number of elements
                insertion_sort(v, l_bound, u_bound);
                return;
            }
            move_pivot_to_last_index(v, l_bound, u_bound);
            int p = partition(v, l_bound, u_bound, u_bound);

            quicksort(v, l_bound, p - 1);
            quicksort(v, p + 1, u_bound);
        }
    }

// Partition (sub-)vector v[l_bound:u_bound] on element with index p
// Returns: index of pivot element after partitioning
    template<typename T>
    inline int partition_fetch_add(std::vector<T> &v, const int l_bound, const int u_bound, const int p, const int block_size,
                                   const int number_of_threads) {
        const int size = u_bound-l_bound+1;
        const int num_blocks = size / block_size;
        const T pivot = v.at(p);
        int return_value;

        if (num_blocks < MINIMUM_BLOCK_NUMBER) {
            return partition(v, l_bound, u_bound, p);
        }

        // indices for synchronized vector access
        std::atomic<int> i(0);
        std::atomic<int> j(0);
        std::atomic<int> k(0);

        std::vector<int> clean_up_left(number_of_threads, -1);
        std::vector<int> clean_up_right(number_of_threads, -1);
        // indices for synchronized accesses to the clean-up vectors
        std::atomic<int> cul(0);
        std::atomic<int> cur(0);

#pragma omp parallel num_threads(number_of_threads) shared(size, num_blocks, pivot, i, j, k)
        {
            bool fetch_left = true;
            bool fetch_right = true;

            int l, r;

            while (std::atomic_fetch_add(&i, 1) < num_blocks) {
                if (fetch_left) {
                    l = std::atomic_fetch_add(&j, 1);
                    fetch_left = false;
                } else if (fetch_right) {
                    r = std::atomic_fetch_add(&k, 1);
                    fetch_right = false;
                }

                if (!fetch_left && !fetch_right) {
                    // block offset
                    int a = 0;
                    int b = 0;
                    // global indices of first block element
                    int v_a = l * block_size;
                    int v_b = size - ((r + 1) * block_size);
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
            /*
             * Clean up preparation:
             *  if fetch_left = true and fetch_right = true -> all blocks were processed completely -> no clean up needed
             *  if fetch_left = false -> left block needs to be processed further
             *  if fetch_right = false -> right block needs ...
            */
            if (!fetch_left) {
                clean_up_left.at(std::atomic_fetch_add(&cul, 1)) = l * block_size;
            }
            if (!fetch_right) {
                clean_up_right.at(std::atomic_fetch_add(&cur, 1)) = size - ((r + 1) * block_size);
            }
        }

        /*
         * Clean up 1st step: (sorting maybe two threads)
         *  Swap elements between remaining left and right blocks
         *  -> everything left of j OR everything right of k is now correctly partitioned
         */

        std::sort(clean_up_left.data(), clean_up_left.data() + cul.load());
        std::sort(clean_up_right.data(), clean_up_right.data() + cur.load(), [](int a, int b) {
            return a > b;
        });

        int a_block = 0;
        int b_block = 0;

        int a = 0;
        int b = 0;

        T buffer;
        bool swap_left = false;
        bool swap_right = false;

        int partition_index;

#pragma omp parallel num_threads(2) shared(a_block, b_block, a, b, buffer, swap_left, swap_right, partition_index)
        {
#pragma omp single nowait
            {
                while (a_block < cul.load() && b_block < cur.load()) {
                    while (a < block_size && b < block_size) {
                        while (a < block_size) {
                            if (v.at(clean_up_left.at(a_block) + a) > pivot) {
                                swap_left = true;
                                break;
                            }
                            a++;
                        }
                        while (b < block_size) {
                            if (v.at(clean_up_right.at(b_block) + b) <= pivot) {
                                swap_right = true;
                                break;
                            }
                            b++;
                        }
                        if (swap_left && swap_right) {
                            buffer = v.at(clean_up_left.at(a_block) + a);
                            v.at(clean_up_left.at(a_block) + a) = v.at(clean_up_right.at(b_block) + b);
                            v.at(clean_up_right.at(b_block) + b) = buffer;
                            swap_left = false;
                            swap_right = false;
                        }
                    }
                    if (a == block_size) {
                        a = 0;
                        a_block++;
                    }
                    if (b == block_size) {
                        b = 0;
                        b_block++;
                    }
                }
            }
#pragma omp single nowait
            {
                /*
                * Clean up 2nd step:
                * partition remainder between L and R
                */
                partition_index = (j.load() * block_size < size - (k * block_size) - 1
                                  ) ? partition_pivot(v, (j.load()) * block_size, size - (k * block_size) - 1, pivot)
                                    : -1;
            }
        }

        if (a_block == cul.load()) {
            int l = (partition_index == -1) ? j.load() * block_size : partition_index + 1;

            /*
             * Cleanup 3rd step:
             * Swap elements > pivot starting from pivot border l
             * with elements < pivot from the remaining blocks
             */
            swap_left = false;
            swap_right = false;
            b = 0;

            while (b_block < cur.load()) {
                while (b < block_size) {
                    if (v.at(clean_up_right.at(b_block) + block_size - 1 - b) <= pivot) {
                        swap_right = true;
                        break;
                    }
                    b++;
                }
                while (l < clean_up_right.at(b_block) + block_size - 1 - b) {
                    if (v.at(l) > pivot) {
                        swap_left = true;
                        break;
                    }
                    l++;
                }
                if ((clean_up_right.at(b_block) + block_size - 1 - b) <= l) {
                    break;
                }
                if (swap_right && swap_left) {
                    buffer = v.at(clean_up_right.at(b_block) + block_size - 1 - b);
                    v.at(clean_up_right.at(b_block) + block_size - 1 - b) = v.at(l);
                    v.at(l) = buffer;
                    swap_left = false;
                    swap_right = false;
                }
                if (b == block_size) {
                    b = 0;
                    b_block++;
                }
            }
            return_value = (v.at(l) <= pivot) ? l : l - 1;
        } else {
            int r = (partition_index == -1) ? size - (k * block_size) - 1 : partition_index;

            /*
             * Cleanup 3rd step:
             * Swap elements > pivot starting from pivot border l
             * with elements < pivot from the remaining blocks
             */
            swap_left = false;
            swap_right = false;
            a = 0;

            while (a_block < cul.load()) {
                while (a < block_size) {
                    if (v.at(clean_up_left.at(a_block) + a) > pivot) {
                        swap_right = true;
                        break;
                    }
                    a++;
                }
                while (r > clean_up_left.at(a_block) + a) {
                    if (v.at(r) <= pivot) {
                        swap_left = true;
                        break;
                    }
                    r--;
                }
                if ((clean_up_left.at(a_block) + a) >= r) {
                    break;
                }
                if (swap_right && swap_left) {
                    buffer = v.at(clean_up_left.at(a_block) + a);
                    v.at(clean_up_left.at(a_block) + a) = v.at(r);
                    v.at(r) = buffer;
                    swap_left = false;
                    swap_right = false;
                }
                if (a == block_size) {
                    a = 0;
                    a_block++;
                }
            }


            return_value = (v.at(r) <= pivot) ? r : r - 1;
        }

        return return_value;
    }

    template<typename T>
    inline T quickselect_parallel(std::vector<T> &v, int l_bound, int u_bound, int k, const int block_size, const int number_of_threads) {
        if (l_bound == u_bound) {
            return v.at(l_bound);
        }
        move_pivot_to_last_index(v, l_bound, u_bound);
        int p = partition_fetch_add(v, l_bound, u_bound, u_bound, block_size, number_of_threads);
        if (p == k) {
            return v.at(p);
        } else if (p < k) {
            return quickselect_parallel(v, p + 1, u_bound, k, block_size, number_of_threads);
        } else {
            return quickselect_parallel(v, l_bound, p - 1, k, block_size, number_of_threads);
        }
    }

    template<typename T>
    void quicksort_parallel(std::vector<T> &v, int l_bound, int u_bound, const int block_size, const int number_of_threads) {
        if (u_bound > l_bound) {
            if (u_bound-l_bound < MINIMUM_VECTOR_ELEMENT_NUMBER) {
                quicksort(v, l_bound, u_bound);
                return;
            }

            move_pivot_to_last_index(v, l_bound, u_bound);
            int p = partition_fetch_add(v, l_bound, u_bound, u_bound, block_size, number_of_threads);

#pragma omp parallel sections
            {
#pragma omp section
                {
                    quicksort_parallel(v, l_bound, p - 1, block_size, number_of_threads/2);
                }
#pragma omp section
                {
                    quicksort_parallel(v, p + 1, u_bound, block_size, number_of_threads/2);
                }
            }
        }
    }

    template<typename T>
    inline T quickselect_parallel_wrapper(std::vector<T> &v, const int k, const int block_size, const int number_of_threads) {
        return quickselect_parallel(v, 0, v.size()-1, k, block_size, number_of_threads);
    }

    template<typename T>
    inline void quicksort_parallel_wrapper(std::vector<T> &v, const int block_size, const int number_of_threads) {
        quicksort_parallel(v, 0, v.size()-1, block_size, number_of_threads);
    }
}