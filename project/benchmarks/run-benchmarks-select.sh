#!/bin/bash
#
#SBATCH --job-name=benchmarks-select
#
#SBATCH --nodes=1
#SBATCH --time=1:00:00

srun g++ -march=native -std=c++11 -fopenmp benchmarks-select.cpp -o benchmarks-select
srun ./benchmarks-select
