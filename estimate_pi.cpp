#include <iomanip>
#include <iostream>
#include <omp.h>
#include <cmath>

using namespace std;

int main()
{
	int num_steps = 400000000;
	double width = 1.0 / double(num_steps);
	double sum = 0.0;
	
	double start_time = omp_get_wtime();
	
#pragma omp parallel num_threads(4)
{
	double partial_sum = 0.0; // Partial sum computed in every thread privately to minimize number of atomic writes to shared variables (sum)
	int num_threads = omp_get_num_threads();
	int block_size = num_steps / num_threads;
	int remainder = num_steps % num_threads;
	int thread_id = omp_get_thread_num();
	int last_block_size = block_size;
	int current_block_size = block_size;
	
	// If steps could not be divided evenly, increment block size and assign remaining steps to last thread
	if (remainder) {
		block_size++;
		last_block_size = num_steps % block_size;
	}
	
	if (thread_id < (num_threads - 1))
		current_block_size = last_block_size;

	for (int i = thread_id * block_size; i < (thread_id * block_size + current_block_size); i++) {
		double x = (i + 0.5) * width;
		partial_sum = partial_sum + (1.0 / (1.0 + x * x));
	}
#pragma omp atomic update
	sum = sum + partial_sum;
}

	double pi = sum  * 4 * width;
	double run_time = omp_get_wtime() - start_time;
	
	cout << "pi with " << num_steps << " steps is " << setprecision(17) << pi << " in " << setprecision(6) << run_time << " seconds\n";
	cout << "Error: " << (M_PI - pi) << "\n";
}
