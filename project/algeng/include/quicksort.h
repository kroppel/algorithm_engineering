#ifndef ALG_ENG_FUNCTIONS
#define ALG_ENG_FUNCTIONS

template <typename T>
int partition_pivot(std::vector<T>& v, int l_bound, int u_bound, T pivot);

template <typename T>
int partition_fetch_add(std::vector<T>& v, int size, int p, int block_size);

template <typename T>
void quicksort(std::vector<T>& v, int l_bound, int u_bound);

template <typename T>
void quicksort_parallel(std::vector<T>& v, int l_bound, int u_bound);

#include "quicksort.cpp"

#endif //ALG_ENG_FUNCTIONS
