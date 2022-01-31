#include <iomanip>
#include <iostream>
#include <omp.h>
#include <mutex>
#include <thread>

using namespace std;

mutex sum_mutex;

void pi_workload(int id, int num_threads, int num_steps, double &sum, double width);

int main() {
  int num_steps = 500000000; // amount of rectangles
  double width = 1.0 / double(num_steps); // width of a rectangle
  double sum = 0.0; // for summing up all heights of rectangles
  double start_time = omp_get_wtime(); // wall clock time in seconds
  int num_threads = 3;
  thread threads [num_threads];
  
  for(int i = 0; i < num_threads; ++i) {
  	threads[i] = thread(pi_workload, i, num_threads, num_steps, ref(sum), width);
  }
  for(int i = 0; i < num_threads; ++i) {
  	threads[i].join();
  }
  
  double pi = sum * 4 * width; // compute pi
  double run_time = omp_get_wtime() - start_time;

  cout << "pi with " << num_steps << " steps is " << setprecision(17)
       << pi << " in " << setprecision(6) << run_time << " seconds\n";
}

void pi_workload(int id, int num_threads, int num_steps, double &sum, double width) {
  double local_sum = 0.0;
  //double x = 0.0;
  
  for (int i = id; i < num_steps; i += num_threads) {
    double x = (i + 0.5) * width; // midpoint
    local_sum = local_sum + (1.0 / (1.0 + x * x)); // add new height of a rectangle
  }
  
  lock_guard<mutex> guard(sum_mutex);
  sum += local_sum;
}
