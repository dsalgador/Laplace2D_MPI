# MPI Implementation to solve the 2-Dimensional Laplace problem using Jacobi iteration method

## Authors
* **Josep Castell**, josepcastellqueralt@gmail.com
* **Dani Salgado**, daniel.salgado@e-campus.uab.cat


# Problem description

The theoretical information about the problem can be found [here](ParallelProgr-MPI.pdf).


## Mathematical problem

![GitHub Logo](http://tutorial.math.lamar.edu/Classes/DE/LaplacesEqn_files/eq0005P.gif) 


## Discretization

![Grid2D](http://basor.fcqb.uasnet.mx/grid.png)


# Tasks

## Task 1
Implement a working version of this algorithm (remember that you can assume that *m* is divisible by *N*). The code must be properly documented. You must include all the necessary information for running your code and indicate the test cases you have used.

## Task 2

Implement a working version of one of these mechanisms (you can make assumptions about the matrix dimensions [m, n] and the number of processes N if needed). You must justify your selection, which means that you must discuss the perceived pros and cons of all methods. The code must be properly documented. You must include all the necessary information for running your code and indicate the test cases you have used.


## Task 3

Make a performance anlysis of your program versions using the given hints and the support of the performance analysis tools available in the lab. You must present an organized explanation of this analysis.



# References

* [1] OpenMPI Slides, UAB.
* [2] MPI tutorial and general information https://computing.llnl.gov/tutorials/mpi/.
* [3] Theoretical video-tutorials of MPI, High Performance Computing by Prof. Matthew Jacob,
Department of Computer Science and Automation, IISC Bangalore. For more details
on NPTEL visit http://nptel.iitm.ac.in. The video-tutorials can be found in
YouTube: https://youtu.be/mzfVimVbguQ, https://youtu.be/mb5wV4AqXso.
* [4] Amdahl’s law (definition of speedup). (2017, July 16). In Wikipedia, The Free Encyclopedia.
Retrieved 22:19, November 25, 2017, from https://en.wikipedia.org/w/index.
php?title=Amdahl%27s_law&oldid=790799480.
* [5] Speedup Ratio and Parallel Efficiency, http://www.bu.edu/tech/support/research/
training-consulting/online-tutorials/matlab-pct/scalability/.
* [6] Measuring parallel scaling performance, https://www.sharcnet.ca/help/index.php/
Measuring_Parallel_Scaling_Performance.


# License
