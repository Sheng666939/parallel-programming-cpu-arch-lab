#!/bin/bash
#PBS -N mpi_ntt
#PBS -l nodes=1:ppn=8
#PBS -j oe

set -euo pipefail

cd "$PBS_O_WORKDIR"

# Adjust THREAD_COUNT, output name, and np according to the experiment.
mpic++ main_v3_hybrid_final.cc -O2 -fopenmp -pthread \
  -DTHREAD_COUNT=2 \
  -DORDINARY_MPI_TASK=0 \
  -DCRT_COLLECT_METHOD=1 \
  -DCRT_INTRA_THREAD=1 \
  -o main_v3_hybrid_t2

mpirun -np 4 ./main_v3_hybrid_t2
