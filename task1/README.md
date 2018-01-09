# Documentation task1

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

* **MPI version** `lapFusion_mpi.c`:
```
mpicc -g -lm -fopenmp -o mpi_lapFusion lapFusion_mpi.c 
```

To execute the programs for example when *N=1024* and the number of iterations is *iter = 100* do:

* **Base version** `lapFusion.c`:
```
./lapFusion 1024 100

```

* **MPI version** `lapFusion_mpi.c`:
```
mpirun -np M mpi_lapFusion 1024 100
```
where M is the number of nodes that one wants to use. To run with 4 nodes the command would be:
```
mpirun -np 4 mpi_lapFusion 1024 100
```

## Test cases

Some code tests can be accessed in earlier commits. If one is interested in those, then have a look at the commit comentaries.
