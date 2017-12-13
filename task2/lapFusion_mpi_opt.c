/*
Basic MPI implementation of the Laplace 2D algorithm using the Jacobi Method.
A grid m x n is assumed and if N is the number of nodes to which distribute work,
it has to be satisfied that m (number of rows) is divisible by N.

As a first approach we assume that n = m, so that we have a grid n x n with N divisible
by n.

  Autors:
      Josep Castell,
      Dani Salgado, 
*/

// Libraries
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

//Definitions
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

/*Initialisation of the grid: internal points set to 0
  and boundary conditions initialised according to the PDF of 
  this assignemnt*/
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


/*Given a matrix * in with nrows rows and ncols columns
  prints it on the console in a representative way*/
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

/*
 This is quite the same as the laplace_step function from the original code
 lapFusion.c. Now my_laplace_step is a function called by each Process to update
 its part of the matrix A (stored in 'my_A'). The part of matrix for each process,
 my_A, has 'nrows' rows and 'ncols' columns. We add also the two parameters 
 'rowstart' and 'rowend' that allow us to decide from which row to wich row my_A is updated
 by the Process that is calling my_laplace_step() function.
 */
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

void copy_A_to_temp(float *A, float * *temp, int nrows, int ncols){
    for(int j = 0; j < nrows; j++){
      for(int i = 0; i < ncols; i++){
        (*temp)[j*ncols +i] = A[j*ncols +i];
      }
    }
}

/*
Commands to run the code:
module load gcc/6.1.0
module load mpe2/mpi-1.10.2/2.4.8
mpicc -g -lm -fopenmp -o mpi_lapFusion lapFusion_mpi_opt.c

mpirun -np N mpi_lapFusion_opt n iter_max k
*/

int main(int argc, char** argv)
{   
  // Initalisation of variables
  double t0, tf; /*Initial and final time counters*/

  int n = 4096; /* Size of the grid n x n */
  int iter_max = 1000; /* Number of iterations */
  int k = 1; /*Number of rows shared (send/recieve) by the processes*/
  float *A, *temp; /* Pointers to grid A and temp */
    
  const float tol = 1.0e-5f; /* Tolerance */
  float error= 1.0f; /* Global error variable */   

  int numtasks, rank, tag = 1,rc; /* Nº of processes, process ID, tag, rc */
  int  my_nrows, my_size; /* Number of rows of my_A, dimension of my_A*/
  float *my_A, *my_temp; /* Portion of A carried by each process*/
  float my_error= 1.0f; /* Error for each process*/
  MPI_Status Stat; /* MPI status variable to control the status*/

  int rowstart, rowend, nrows; /*Auxiliar variables related to rows*/

  //INIT MPI environment
  rc = MPI_Init (&argc, &argv);
  if (rc != MPI_SUCCESS)
    {
      printf ("Error starting MPI program. Terminating.\n");
      MPI_Abort (MPI_COMM_WORLD, rc);
      return -1;
    }
  MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank); 
  //END basic INIT MPI environment


  // Abort the program if the number of processes is less than 2
  if(numtasks < 2){
    printf ("This program works with 2 or more processes (-np N with N >=2).\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
    return -1;
  }


  //BEGIN MASTER Initialisation of A, temp and initial time
  if(rank == MASTER){
  t0 = MPI_Wtime(); //Record the initial time 

  // get runtime arguments 
  if (argc>1) {  n        = atoi(argv[1]); }
  if (argc>2) {  iter_max = atoi(argv[2]); }
  if (argc>3) {  k = atoi(argv[3]); }

  // Allocate memory for A and temp
  if( ( A = (float*) malloc(n*n*sizeof(float)) ) == NULL ){
    printf ("Error when allocating memory for A.\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
    return -1;
  }
  if( ( temp = (float*) malloc(n*n*sizeof(float)) ) == NULL ){
    printf ("Error when allocating memory for temp.\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
    return -1;
  }

  // set boundary conditions
  laplace_init (A, n);
  laplace_init (temp, n);
  A[(n/128)*n+n/128] = 1.0f; // set singular point

  printf("Jacobi relaxation Calculation: %d x %d mesh,"
         " maximum of %d iterations\n", 
         n, n, iter_max );
  } //END MASTER initialisation

  //All processes initialise iter to 0
  int iter = 0;

  //Broadcast de global (MASTER) error, n, iter_max and k
  MPI_Bcast(&error, 1, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&n, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&iter_max, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&k, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

  //Initialise some auxiliar variables
  my_nrows = n/numtasks; 
  nrows = my_nrows +2*k;
  my_size = n*nrows;

  //Allocate memory for my_A and my_temp
  if( ( my_A = (float*) malloc( my_size*sizeof(float)) ) == NULL ){
    printf ("Error when allocating memory for my_A.\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
    return -1;
  }
  if( ( my_temp = (float*) malloc(my_size*sizeof(float)) ) == NULL ){
    printf ("Error when allocating memory for my_temp.\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
    return -1;
  }

  int starting = n*k;
  //int ending = n*k;
  //Distribute the rows of A and temp among all the processes --> store in my_A, my_temp
  MPI_Scatter(A, my_nrows*n,  MPI_FLOAT, my_A+starting, my_nrows*n, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
  //float * my_temp_plusnk = my_temp+starting;
  //copy_A_to_temp(my_A+starting, &my_temp_plusnk, my_nrows, n);
  MPI_Scatter(temp, my_nrows*n,  MPI_FLOAT, my_temp+starting, my_nrows*n, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
 

 while ( error > tol*tol && iter < iter_max )
  {
    iter++;
    /*Send and Recv calls so that each process obtain two additional rows needed for
    the computation of the new values. 
    */   
    
    if( (iter % k) == 0 )
    {    
      if(rank > MASTER){ 
        //For all the processes apart from MASTER, which does not need a previous row

        /* Send the first k rows of the process 'rank' to the process 'rank-1'*/
        MPI_Send(my_A+n*k, n*k, MPI_FLOAT, rank-1, tag ,MPI_COMM_WORLD);
        /* Process 'rank' recieves the last k rows from the process 'rank-1'*/
        MPI_Recv(my_A  , n*k, MPI_FLOAT, rank-1, tag ,MPI_COMM_WORLD, &Stat);
      }
      if(rank < numtasks -1 ){
         //For all the processes apart from THE LAST, which does not need a 'last' row

         /* Process 'rank' recieves the first k rows from the process 'rank+1'*/
         MPI_Recv( (my_A + n*(my_nrows-k+2) )  , n*k, MPI_FLOAT, rank+1, tag ,MPI_COMM_WORLD, &Stat);
         /* Send the last k rows of the process 'rank' to the process 'rank+1'*/
         MPI_Send(  (my_A + n*(my_nrows-k+1))  , n*k, MPI_FLOAT, rank+1, tag ,MPI_COMM_WORLD);       
      }

    } //endif iter % k
    
    /*Set values for rowstart and rowend in order to make all processes modify only the internal
    points of A. We have to distinguish between the process that has the first block of rows of A
    (the MASTER) and the one that has the last block of rows (the process numtasks-1).
    */
    if(rank == MASTER){
      rowstart = 1+k;//2
      rowend = nrows-1;//nrows -1 
    }
    else if(rank == (numtasks - 1)){
      rowstart = k; //1
      rowend = nrows - (1+k); // nrows -2 
    }
    else{
      rowstart = k; //1
      rowend = nrows-k; //nrows-1   
    }
    //Each process performs the laplace_step updating the points of my_A that are interior points of A
    my_error= my_laplace_step(my_A, my_temp, nrows, n, rowstart, rowend);
    
    /*Reduction operation: the maximum among all my_error from all processes is calculated and stored
    in the variable error, which is the global error and originally stored in the MASTER process*/
    MPI_Reduce(&my_error, &error, 1, MPI_FLOAT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    //Swap the roles of my_A and my_temp (double buffer) to be prepared for the next iteration.
    float *swap= my_A; my_A=my_temp; my_temp= swap;
  }
  /*The master process gathers all the final portions of A stored in my_A of each process to build
   the matrix A corresponding to the last iteration.
  */
  //MPI_Scatter(A, my_nrows*n,  MPI_FLOAT, my_A+starting, my_nrows*n, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
 
  MPI_Gather(my_A+starting, my_nrows*n ,  MPI_FLOAT, A, my_nrows*n, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
  
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
