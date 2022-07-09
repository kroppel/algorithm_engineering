#ifndef ALGENG
#define ALGENG
#include <omp.h>
#include <iostream>
#include <vector>
#include <atomic>
#include <functional>
#include <iterator>
#include <algorithm>

namespace algeng {
    template<typename T>
    inline void median_of_3(std::vector<T> &v, int l_bound, int u_bound) {
        // Compute median of 3 elements that are distributed evenly over vector and swap
        // the element with the element on index 'u_bound'.
        // Worst Case O(N^2) Search for Median, but it is only 3 elements.
        int p1[3] = {l_bound, (u_bound+l_bound)/2, u_bound};
        for (int a:p1) {
            int sorted_index=0;
            for (int b:p1) {
                if (v.at(a)>v.at(b))
                    sorted_index++;
            }
            if (sorted_index==1) {
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

    template<typename T>
    inline void median_of_5(std::vector<T> &v, int l_bound, int u_bound) {
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

// Partition (sub-)vector v[l_bound:u_bound] on pivot element on index u_bound
// Returns: index of pivot element after partitioning
    template<typename T>
    inline int partition(std::vector<T> &v, int l_bound, int u_bound) {
        T buffer;
        int i = l_bound;
        int j = u_bound - 1;

        if (u_bound > l_bound) {
            T pivot = v.at(u_bound);

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
                v.at(u_bound) = v.at(i);
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
            return i-1;
        }
        else if (u_bound == l_bound) {
            if (v.at(u_bound) <= pivot)
                return u_bound;
            else
                return u_bound-1;
        }
        else {
            return u_bound;
        }
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
    inline void insertion_sort_descending(std::vector<T> &v, int l_bound, int u_bound) {
        int j;
        T current;
        for (int i = l_bound; i <= u_bound; i++) { // i=3
            current = v.at(i); //1
            j = i - 1; //j=-1

            while (j <= l_bound && v.at(j) > current) {
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
            median_of_3(v, l_bound, u_bound);
            int p = partition(v, l_bound, u_bound);

            quicksort(v, l_bound, p - 1);
            quicksort(v, p + 1, u_bound);
        }
    }

// Partition (sub-)vector v[l_bound:u_bound] on element with index p
// Returns: index of pivot element after partitioning
    
    /**
     * Partition (partial) vector multithreaded.
     *
     * Partition (partial) vector using the specified amound of threads and block size.
     *
     * @param &v Reference of input vector.
     * @param l_bound Lower bound defining the part of the vector to be partitioned.
     * @param u_bound Upper bound defining the part of the vector to be partitioned.
     * @param block_size Block size to be used for partitioning.
     * @param number_of_threads Number of threads to be used for partitioning.
     */
    template<typename T>
    inline int partition_fetch_add(std::vector<T> &v, const int l_bound, const int u_bound, const int block_size,
                                   const int number_of_threads) {
        // size leaves out pivot element at end of vector because it will not be included in the processing
        const int size = u_bound-l_bound;
        const int num_blocks = size / block_size;
        const T pivot = v.at(u_bound);
        int return_value;

        if (num_blocks < number_of_threads) {
            return partition(v, l_bound, u_bound);
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

#pragma omp parallel num_threads(number_of_threads) shared(i, j, k)
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
                    int v_a = l_bound + l * block_size;
                    int v_b = l_bound + size - ((r + 1) * block_size);
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
                clean_up_left.at(std::atomic_fetch_add(&cul, 1)) = l_bound + l * block_size;
            }
            if (!fetch_right) {
                clean_up_right.at(std::atomic_fetch_add(&cur, 1)) = l_bound + size - ((r + 1) * block_size);
            }
        }

        /*
         * Clean up 1st step: (sorting maybe two threads)
         *  Swap elements between remaining left and right blocks
         *  -> everything left of j OR everything right of k is now correctly partitioned
         */

        /*insertion_sort(clean_up_left, 0, cul.load());
        insertion_sort_descending(clean_up_right, 0, cur.load());*/

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
                partition_index = (l_bound + j.load() * block_size <= l_bound + size - (k.load() * block_size) - 1
                                  ) ? partition_pivot(v, l_bound + (j.load()) * block_size, l_bound + size - (k.load() * block_size) - 1, pivot)
                                    : -1;
            }
        }

        if (a_block == cul.load()) {
            int l = (partition_index == -1) ? l_bound + j.load() * block_size : partition_index + 1;

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

            return_value = (v.at(l) > pivot) ? l: l+1;
            buffer = v.at(u_bound);
            v.at(u_bound) = v.at(return_value);
            v.at(return_value) = buffer;

        } else {
            int r = (partition_index == -1) ? l_bound + size - (k.load() * block_size) - 1 : partition_index;

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
                        swap_left = true;
                        break;
                    }
                    a++;
                }
                while (r > clean_up_left.at(a_block) + a) {
                    if (v.at(r) <= pivot) {
                        swap_right = true;
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

            return_value = (v.at(r) > pivot) ? r: r+1;
            buffer = v.at(u_bound);
            v.at(u_bound) = v.at(return_value);
            v.at(return_value) = buffer;
        }

        return return_value;
    }
    
    /**
     * Select kth vector element multithreaded.
     *
     * Select the kth element from (partial) vector multithreaded using the specified number of threads
     * and block size.
     *
     * @param &v Reference of vector to be sorted.
     * @param l_bound Lower bound defining the part of the vector to be selected from.
     * @param u_bound Upper bound defining the part of the vector to be selected from.
     * @param k Position of the element to be retrieved in vector v if v was sorted.
     * @param block_size Block size to be used for element selection.
     * @param number_of_threads Number of threads to be used for element selection.
     */
    template<typename T>
    inline T quickselect_parallel(std::vector<T> &v, int l_bound, int u_bound, int k, const int block_size, const int number_of_threads) {
        if (l_bound == u_bound) {
            return v.at(l_bound);
        }
        median_of_3(v, l_bound, u_bound);
        int p = partition_fetch_add(v, l_bound, u_bound, block_size, number_of_threads);
        if (p == k) {
            return v.at(p);
        } else if (p < k) {
            return quickselect_parallel(v, p + 1, u_bound, k, block_size, number_of_threads);
        } else {
            return quickselect_parallel(v, l_bound, p - 1, k, block_size, number_of_threads);
        }
    }

    /**
     * Sort vector multithreaded.
     *
     * Sort the given partial vector multithreaded using the specified number of threads
     * and block size.
     *
     * @param &v Reference of vector to be sorted.
     * @param l_bound Lower bound defining the part of the vector to be sorted.
     * @param u_bound Upper bound defining the part of the vector to be sorted.
     * @param block_size Block size to be used for sorting.
     * @param number_of_threads Number of threads to be used for sorting.
     */
    template<typename T>
    void quicksort_parallel(std::vector<T> &v, int l_bound, int u_bound, const int block_size, const int number_of_threads) {
        if (u_bound > l_bound) {
            if (u_bound-l_bound < 100000 or number_of_threads<=1) {
                
                quicksort(v, l_bound, u_bound);
                return;
            }

            else if (u_bound-l_bound < 1000000) {
#pragma omp parallel 
                {
#pragma omp single nowait
                    {
                        median_of_3(v, l_bound, u_bound);
                        int p = partition(v, l_bound, u_bound);
#pragma omp task
                        {
                            quicksort_parallel(v, l_bound, p, block_size, number_of_threads*(p-l_bound)/(u_bound-l_bound));
                        }
#pragma omp task
                        {
                            quicksort_parallel(v, p + 1, u_bound, block_size, number_of_threads*(u_bound-p)/(u_bound-l_bound));
                        }
                    }
#pragma omp taskwait
                }
            }

            else {
                median_of_3(v, l_bound, u_bound);
                int p = partition_fetch_add(v, l_bound, u_bound, block_size, number_of_threads);

                quicksort_parallel(v, l_bound, p-1, block_size, number_of_threads);
                quicksort_parallel(v, p + 1, u_bound, block_size, number_of_threads);
            }
        }
    }

    /**
     * Select kth vector element multithreaded.
     *
     * Select the kth vector element multithreaded using the specified number of threads
     * and block size.
     *
     * @param &v Reference of vector to be selected from.
     * @param k Position of the element to be retrieved in vector v if v was sorted.
     * @param block_size Block size to be used for element selection.
     * @param number_of_threads Number of threads to be used for element selection.
     */
    template<typename T>
    inline T quickselect_parallel_wrapper(std::vector<T> &v, const int k, const int block_size, const int number_of_threads) {
        return quickselect_parallel(v, 0, v.size()-1, k, block_size, number_of_threads);
    }

    /**
     * Sort vector multithreaded.
     *
     * Sort the given vector multithreaded using the specified number of threads
     * and block size.
     *
     * @param &v Reference of vector to be sorted.
     * @param block_size Block size to be used for sorting.
     * @param number_of_threads Number of threads to be used for sorting.
     */
    template<typename T>
    inline void quicksort_parallel_wrapper(std::vector<T> &v, const int block_size, const int number_of_threads) {
        quicksort_parallel(v, 0, v.size()-1, block_size, number_of_threads);
    }
}

#endif //ALGENG
