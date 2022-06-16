#include "quicksort.hpp"
#include <random>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <parallel/algorithm>
#include <string>

using namespace algeng;

int main() {
    double start_time, runtime1, runtime2, runtime3;

    const size_t sizes[4] = {1000000, 10000000, 20000000, 50000000};

    const int number_threads = 48;

    for (const size_t size:sizes) {

        const int block_size = (int) size*0.00005;

        std::ofstream outfile;
        outfile.open("/home/hk-project-test-scalfsu/hgf_qei8127/code/benchmark-select-" + std::to_string(size) + ".txt");

        std::vector<int> v1, v2, v3;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> d(0, size);
        // Fill vectors with random numbers between 0 and size - 1
        outfile << "Vector size " << size << "\n";
        outfile << "Block size " << block_size << "\n";
        outfile << "Number Threads " << number_threads << "\n\n";
        for (int i = 0; i < size; ++i) {
            v1.push_back(d(gen));
        }

        v2 = v1;
        v3 = v1;

        int m = size/2;

        start_time = omp_get_wtime();
        std::nth_element(v1.begin(), v1.begin() + m, v1.end());
        runtime1 = omp_get_wtime() - start_time;

        start_time = omp_get_wtime();
        int element_n_test = quickselect_parallel_wrapper(v2, m, block_size, number_threads);
        runtime2 = omp_get_wtime() - start_time;

        start_time = omp_get_wtime();
        __gnu_parallel::nth_element(v3.begin(), v3.begin() + m, v3.end());
        runtime3 = omp_get_wtime() - start_time;

        outfile << "Is nth element: " << (v1.at(m) == element_n_test) << "\n";

        outfile << "std::nth_element with runtime: " << std::setprecision(6) << runtime1 << "s\n";
        outfile << "quickselect_parallel_wrapper with runtime: " << std::setprecision(6) << runtime2 << "s\n";
        outfile << "__gnu_parallel::nth_element with runtime: " << std::setprecision(6) << runtime3 << "s\n";

        outfile.close();
    }
}
