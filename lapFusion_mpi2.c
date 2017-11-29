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

float laplace_step(float *in, float *out, int n)
{
  int i, j;
  float error=0.0f;
  for ( j=1; j < n-1; j++ )
    #pragma omp simd reduction(max:error)
    for ( i=1; i < n-1; i++ )
    {
      out[j*n+i]= stencil(in[j*n+i+1], in[j*n+i-1], in[(j-1)*n+i], in[(j+1)*n+i]);
      error = max_error( error, out[j*n+i], in[j*n+i] );
    }
  return error;
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

int main(int argc, char** argv)
{
  int n = 4096;
  int iter_max = 1000;
  float *A, *temp;
    
  const float tol = 1.0e-5f;
  float error= 1.0f;   

  int numtasks, rank, dest, source, tag = 1,rc;
  int ri, rf, my_nrows, posi, posf, my_size;
  float *my_A, *my_temp;
  float my_error= 1.0f;
  MPI_Status Stat;


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

  if(rank == MASTER){

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
  
  if(n<20){
    print_matrix(A, n,n); 
 }

  //send the parts of matrix to each process

  } 

  int iter = 0;

  if(numtasks >1){
  //Broadcast de global error, n and iter_max
  MPI_Bcast(&error, 1, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&n, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&iter_max, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  //

  my_nrows = n/numtasks; 
  my_size = n*(my_nrows+2);

  my_A    = (float*) malloc( my_size*sizeof(float) );
  my_temp = (float*) malloc( my_size*sizeof(float) );

  MPI_Scatter(A, my_nrows*n,  MPI_FLOAT, my_A+n, my_nrows*n, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
  MPI_Scatter(temp, my_nrows*n,  MPI_FLOAT, my_temp+n, my_nrows*n, MPI_FLOAT, MASTER, MPI_COMM_WORLD);


  ri = rank * my_nrows; //-1
  rf = ri + my_nrows;//+1
  posi = ri*n;
  posf = rf*n;



 }


  /*
  while ( error > tol*tol && iter < iter_max )
  {
    iter++;
    error= laplace_step (A, temp, n);
    float *swap= A; A=temp; temp= swap; // swap pointers A & temp
  }
*/
  if(rank == MASTER){
    error = sqrtf( error );
    printf("Total Iterations: %5d, ERROR: %0.6f, ", iter, error);
    printf("A[%d][%d]= %0.6f\n", n/128, n/128, A[(n/128)*n+n/128]);

    free(A); free(temp);
   }
   MPI_Finalize();

   return 0;

}
