#!/bin/bash
module load openmpi/1.8.1
module load gcc/6.1.0

N=$1
n=$2
iter=$3
k=$4

mpirun -np $N /home/master/ppM/ppM-1-7/Escritorio/EntregaMPI/task3/mpi_lapFusion_opt $n $iter $k #>> StrongScalingSingleNode.csv;
