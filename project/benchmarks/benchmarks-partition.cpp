#include "quicksort.hpp"
#include <random>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <parallel/algorithm>
#include <string>

using namespace algeng;

void benchmark_one_iteration() {
    double start_time, runtime1, runtime2, runtime3;

    const size_t sizes[5] = {1000000, 10000000, 100000000, 200000000, 500000000};

    const int number_threads = 48;

    for (const size_t size:sizes) {
        const int block_size = (int) size*0.00005;

        std::ofstream outfile;
        outfile.open("/home/hk-project-test-scalfsu/hgf_qei8127/code/benchmark-partition-" + std::to_string(size) + ".txt");

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

        int pivot = v1.at(size-1);

        start_time = omp_get_wtime();
        std::partition(v1.begin(), v1.end(), [&pivot](int x) {return x<=pivot;});
        runtime1 = omp_get_wtime() - start_time;

        start_time = omp_get_wtime();
        partition_fetch_add(v2, 0, v2.size()-1, block_size, number_threads);
        runtime2 = omp_get_wtime() - start_time;

        start_time = omp_get_wtime();
        __gnu_parallel::partition(v3.begin(), v3.end(), [&pivot](int x) {return x<=pivot;});
        runtime3 = omp_get_wtime() - start_time;

        outfile << "Is partitioned: " << std::is_partitioned(v2.begin(), v2.end(), [&pivot](int x) {return x<=pivot;}) << "\n";

        outfile << "std::partition with runtime: " << std::setprecision(6) << runtime1 << "s\n";
        outfile << "partition_fetch_add with runtime: " << std::setprecision(6) << runtime2 << "s\n";
        outfile << "__gnu_parallel::partition with runtime: " << std::setprecision(6) << runtime3 << "s\n";

        outfile.close();
    }
}

int main() {
    benchmark_one_iteration();
}
