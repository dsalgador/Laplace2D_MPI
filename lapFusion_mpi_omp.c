#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MASTER 0        /* task ID of master task */

float stencil ( float v1, float v2, float v3, float v4)
{
  return (v1 + v2 + v3 + v4) * 0.25f;
}

float max_error ( float prev_error, float old, float new )
{
  float t= fabsf( new - old );
  return t>prev_error? t: prev_error;
}

void laplace_init(float *in, int n)
{
  int i;
  const float pi  = 2.0f * asinf(1.0f);
  memset(in, 0, n*n*sizeof(float));
  for (i=0; i<n; i++) {
    float V = in[i*n] = sinf(pi*i / (n-1));
    in[ i*n+n-1 ] = V*expf(-pi);
  }
}

void print_matrix(float * in, int nrows, int ncols){
    int i,j;
    for ( j=0; j < nrows; j++ ){
      for ( i=0; i < ncols; i++ )
    {
      printf("%0.3f|", in[j*ncols+i]);
    }
    printf("\n");
    }
    printf("\n");
}

float my_laplace_step(float *in, float *out, int nrows, int ncols, int rowstart, int rowend)
{
  int i, j;
  float my_error=0.0f;
  for ( j=rowstart; j < rowend; j++ )
    #pragma omp simd reduction(max:my_error)
    for ( i=1; i < ncols-1; i++ )
    {
      out[j*ncols+i]= stencil(in[j*ncols+i+1], in[j*ncols+i-1], in[(j-1)*ncols+i], in[(j+1)*ncols+i]);
      my_error = max_error( my_error, out[j*ncols+i], in[j*ncols+i] );
    }
  return my_error;
}

/*
Commands to run the code:
module load gcc/6.1.0
module load mpe2/mpi-1.10.2/2.4.8
mpicc -g -lm -fopenmp -o mpi_lapFusion2 lapFusion_mpi2.c

mpirun -np N mpi_lapFusion2 n iter_max
*/

int main(int argc, char** argv)
{  
  double t0, tf;

  int n = 4096;
  int iter_max = 1000;
  float *A, *temp;
    
  const float tol = 1.0e-5f;
  float error= 1.0f;   

  int numtasks, rank, tag = 1,rc;
  int  my_nrows, my_size;
  float *my_A, *my_temp;
  float my_error= 1.0f;
  MPI_Status Stat;

  int rowstart, rowend, nrows;

  //INIT MPI ENVIRONMENT
  rc = MPI_Init (&argc, &argv);
  if (rc != MPI_SUCCESS)
    {
      printf ("Error starting MPI program. Terminating.\n");
      MPI_Abort (MPI_COMM_WORLD, rc);
      return -1;
    }
  MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank); 
  
  //Abort the program if the number of processes is less than 2
  if(numtasks < 2){
    printf ("This program works with 2 or more processes (-np N with N >=2).\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
    return -1;
  }


  //Initialisation of A, temp and initial time
  if(rank == MASTER){
  t0 = MPI_Wtime(); //Record the initial time 

  // get runtime arguments 
  if (argc>1) {  n        = atoi(argv[1]); }
  if (argc>2) {  iter_max = atoi(argv[2]); }

  A    = (float*) malloc( n*n*sizeof(float) );
  temp = (float*) malloc( n*n*sizeof(float) );

  //  set boundary conditions
  laplace_init (A, n);
  laplace_init (temp, n);
  A[(n/128)*n+n/128] = 1.0f; // set singular point

  printf("Jacobi relaxation Calculation: %d x %d mesh,"
         " maximum of %d iterations\n", 
         n, n, iter_max );
  } 

  //All processes initialise iter to 0
  int iter = 0;

  //Broadcast de global (MASTER) error, n and iter_max
  MPI_Bcast(&error, 1, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&n, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&iter_max, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

  //Initialise some variables
  my_nrows = n/numtasks; 
  nrows = my_nrows +2;
  my_size = n*(my_nrows+2);

  //Alloc memory for my_A and my_temp
  my_A    = (float*) malloc( my_size*sizeof(float) );
  my_temp = (float*) malloc( my_size*sizeof(float) );

  //Distribute the rows of A and temp among all the processes --> store in my_A, my_temp
  MPI_Scatter(A, my_nrows*n,  MPI_FLOAT, my_A+n, my_nrows*n, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
  MPI_Scatter(temp, my_nrows*n,  MPI_FLOAT, my_temp+n, my_nrows*n, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
 

 while ( error > tol*tol && iter < iter_max )
  {
    iter++;
    /*Send and Recv calls so that each process obtain two additional rows needed for
    the computation of the new values. 
    */
    if(rank > MASTER){
      MPI_Send(my_A+n, n, MPI_FLOAT, rank-1, tag ,MPI_COMM_WORLD);
      MPI_Recv(my_A  , n, MPI_FLOAT, rank-1, tag ,MPI_COMM_WORLD, &Stat);
    }
    if(rank < numtasks -1 ){
       MPI_Send(  (my_A + n*(my_nrows))  , n, MPI_FLOAT, rank+1, tag ,MPI_COMM_WORLD);
       MPI_Recv( (my_A + n*(my_nrows+1) )  , n, MPI_FLOAT, rank+1, tag ,MPI_COMM_WORLD, &Stat);
    }

    /*Set values for rowstart and rowend in order to make all processes modify only the internal
    points of A. We have to distinguish between the process that has the first block of rows of A
    (the MASTER) and the one that has the last block of rows (the process numtasks-1).
    */
    if(rank == MASTER){
      rowstart =2;
      rowend = nrows-1;
    }
    else if(rank == (numtasks - 1)){
      rowstart = 1;
      rowend = nrows -2;
    }
    else{
      rowstart = 1;
      rowend = nrows-1;   
    }
    //Each process perform the laplace_step updating the points of my_A that are interior points of A
    my_error= my_laplace_step(my_A, my_temp, nrows, n, rowstart, rowend);
    
    /*Reduction operation: the maximum among all my_error from all processes is calculated and stored
    in the variable error, which is the global error and originally stored in the MASTER process*/
    MPI_Reduce(&my_error, &error, 1, MPI_FLOAT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    //Swap the roles of my_A and my_temp (double buffer) to be prepared for the next iteration.
    float *swap= my_A; my_A=my_temp; my_temp= swap;
  }
  /*The master process gather all the final portions of A stored in my_A of each process to build
   the matrix A corresponding to the last iteration.
  */
  MPI_Gather(my_A+n, my_nrows*n ,  MPI_FLOAT, A, my_nrows*n, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
  
  /*The MASTER process computes the final error as the sqrt of the variable error
    and some information is printed onto the screen*/
  if(rank == MASTER){
    error = sqrtf( error );
    printf("Total Iterations: %5d, ERROR: %0.6f, ", iter, error);
    printf("A[%d][%d]= %0.6f\n", n/128, n/128, A[(n/128)*n+n/128]);
    free(A); free(temp);
   }
   /*The master process prints the execution time*/
   if(rank == MASTER){
    tf = MPI_Wtime();
    printf("Elapsed time, %2.5lf\n", tf-t0);
   }
   //Finalize the MPI environment.
   MPI_Finalize();
   return 0;
}
