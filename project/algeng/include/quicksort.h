#ifndef ALG_ENG_FUNCTIONS
#define ALG_ENG_FUNCTIONS

template <typename T>
void quicksort(std::vector<T>& v, int l_bound, int u_bound);

template <typename T>
void quicksort_parallelized(std::vector<T>& v, int l_bound, int u_bound);

#include "quicksort.cpp"

#endif //ALG_ENG_FUNCTIONS
