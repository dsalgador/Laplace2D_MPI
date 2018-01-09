# Documentation task2

## Prerequisites
Some external libraries will need to be included to run this program. 
### Libraries, modules and compiler

The MPI version program, `lapFusion_mpi.c` supports the following C compiler versions:

* Linux gcc v 5.4.0 or higher
* Linux openmpi v 1.8.1 (or higher)

The base version program `lapFusion.c` only needs some of the gcc versions above.

Moreover we need to include the following modules:

```
module load gcc/6.1.0
module load openmpi/1.8.1
```

### Compilation and execution

To compile each code version

* **Base version** `lapFusion.c`:
```
gcc -fopenmp -lm -Wall lapFusion.c -o lapFusion

```

* **Optimized MPI version** `lapFusion_mpi_opt.c`:
```
mpicc -g -lm -fopenmp -o mpi_lapFusion lapFusion_mpi_opt.c 
```

To execute the program for example when *N=1024* and when the number of iterations is *iter = 100* do:

* **Base version** `lapFusion.c`:
```
./lapFusion 1024 100

```

* **Optimized MPI version** `lapFusion_mpi_opt.c`:
To execute the program for example when *N=1024*, the number of iterations is *iter = 100* and the number of shared rows is
*k = 4* (or iterations between comunications) do:
```
mpirun -np M mpi_lapFusion 1024 100 4
```
where M is the number of processes that one wants to use. To run with 2 processes the command would be:
```
mpirun -np 2 mpi_lapFusion 1024 100 4
```