#include <algorithm>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <random>
using namespace std;

vector<int> get_random_int_vector(int n) {
  default_random_engine re{random_device{}()};
  uniform_int_distribution<int> next_rand{INT32_MIN, INT32_MAX};
  vector<int> v(n);
  for (auto &num : v) {
    num = next_rand(re);
  }
  return v;
}

void insertion_sort(int *arr, int n) {
    for (int i = 1; i < n; i++) {
        int j = i;
        int buff = arr[j];
        while ((j > 0) && (buff < arr[j - 1])) {
            arr[j] = arr[j - 1];
            j--;
        }
        arr[j] = buff;
    }
}

inline void merge(const int *__restrict__ a, const int *__restrict__ b,
                  int *__restrict__ c, const int a_size, const int b_size,
                  const int c_size) {
  int idx_a = 0;
  int idx_b = 0;
  for (int i = 0; i < c_size; ++i) {
    if (idx_a == a_size) {
      c[i] = b[idx_b++];
    } else if (idx_b == b_size) {
      c[i] = a[idx_a++];
    } else {
      c[i] = (a[idx_a] < b[idx_b]) ? a[idx_a++] : b[idx_b++];
    }
  }
}

void merge_sort_speed_up(int *arr, int n) {
    // use insertion sort for small n
    if (n < 100) {
        insertion_sort(arr, n);
    }
    else {
        const int size_a = n / 2;
        const int size_b = n - size_a;
        // make next recursive call a task
#pragma omp task
        merge_sort_speed_up(arr, size_a); // recursive call
        merge_sort_speed_up(arr + size_a, size_b); // recursive call
        // here should be a taskwait
#pragma omp taskwait
        if (n > 400) {
            int *c = new int[n];
            merge(arr, arr + size_a, c, size_a, size_b, n);
            memcpy(arr, c, sizeof(int) * n);
            delete[](c);
        }
        else {
            int c[n];
            merge(arr, arr + size_a, c, size_a, size_b, n);
            memcpy(arr, c, sizeof(int) * n);
        }
    }
}

void merge_sort_parallel(int *arr, int n) {
#pragma omp parallel
#pragma omp single nowait
    merge_sort_speed_up(arr, n);
}

int main(int argc, char *argv[]) {
  const int n = 10000000;
  vector<int> v = get_random_int_vector(n);
  vector<int> v_copy = v;

  double start = omp_get_wtime();
  merge_sort_parallel(v.data(), n);
  cout << "speed_up: " << omp_get_wtime() - start << " seconds" << endl;

  start = omp_get_wtime();
  sort(begin(v_copy), end(v_copy));
  cout << "std::sort: " << omp_get_wtime() - start << " seconds" << endl;

  if (v != v_copy) {
    cout << "sort implementation is buggy\n";
  }
}
