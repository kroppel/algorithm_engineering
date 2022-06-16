#!/bin/bash
#
#SBATCH --job-name=benchmarks-sort
#
#SBATCH --nodes=1
#SBATCH --time=1:00:00

srun g++ -march=native -std=c++11 -fopenmp benchmarks-partition.cpp -o benchmarks-partition
srun ./benchmarks-partition
