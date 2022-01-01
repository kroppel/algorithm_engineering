#include <iomanip>
#include <iostream>
#include <omp.h>
#include <mutex>

using namespace std;

int main() {
  mutex sum_mutex;
  int num_steps = 100000000; // amount of rectangles
  double width = 1.0 / double(num_steps); // width of a rectangle
  double start_time = omp_get_wtime(); // wall clock time in seconds
  double pi = 0.0;
#pragma omp parallel num_threads(4)
  {
    int num_threads = omp_get_num_threads();
    int thread_id = omp_get_thread_num();
    double sum_local = 0.0; // for summing up heights locally
    for (int i = thread_id; i < num_steps; i += num_threads) {
      double x = (i + 0.5) * width; // midpoint
      sum_local = sum_local + (1.0 / (1.0 + x * x)); // add new height
    }
    //sum_mutex.lock();
    pi += sum_local * 4 * width;
    //sum_mutex.unlock();
  }
  double run_time = omp_get_wtime() - start_time;
  cout << "pi with " << num_steps << " steps is " << setprecision(17) 
	    << pi << " in " << setprecision(6) << run_time << " seconds\n";
}
