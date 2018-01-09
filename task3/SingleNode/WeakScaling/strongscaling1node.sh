#!/bin/bash
#module load openmpi/1.8.1
#module load gcc/6.1.0

N=$1
n=$2
iter=$3
k=$4
filename1=$5
#filename2=$6

if [ $N == 1 ] 
then
	/home/master/ppM/ppM-1-7/Escritorio/EntregaMPI/task3/SingleNode/lapFusion $n $iter >> $filename1 #
else
	mpirun -np $N /home/master/ppM/ppM-1-7/Escritorio/EntregaMPI/task3/mpi_lapFusion_opt $n $iter $k  >> $filename1 #>> StrongScalingSingleNode.csv;
fi