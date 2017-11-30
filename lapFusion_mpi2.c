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

float my_laplace_step(float *in, float *out, int nrows, int ncols)
{
  int i, j;
  float my_error=0.0f;
  for ( j=1; j < nrows-1; j++ )
    #pragma omp simd reduction(max:my_error)
    for ( i=1; i < ncols-1; i++ )
    {
      out[j*ncols+i]= stencil(in[j*ncols+i+1], in[j*ncols+i-1], in[(j-1)*ncols+i], in[(j+1)*ncols+i]);
      my_error = max_error( my_error, out[j*ncols+i], in[j*ncols+i] );
    }
  return my_error;
}

/*
Commands to test this version
module load gcc/6.1.0
module load mpe2/mpi-1.10.2/2.4.8
mpicc -g -lm -o mpi_lapFusion2 lapFusion_mpi2.c

mpirun -np 3 mpi_lapFusion2 12 2
*/

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

  //if(numtasks >1){
  //Broadcast de global error, n and iter_max
  MPI_Bcast(&error, 1, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&n, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&iter_max, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  //

  my_nrows = n/numtasks; 
  my_size = n*(my_nrows+2);

  if (rank == MASTER){
    printf("my_nrows %d, process %d\n", my_nrows ,rank);
  }


  my_A    = (float*) malloc( my_size*sizeof(float) );
  my_temp = (float*) malloc( my_size*sizeof(float) );

  MPI_Scatter(A, my_nrows*n,  MPI_FLOAT, my_A+n, my_nrows*n, MPI_FLOAT, MASTER, MPI_COMM_WORLD);
  MPI_Scatter(temp, my_nrows*n,  MPI_FLOAT, my_temp+n, my_nrows*n, MPI_FLOAT, MASTER, MPI_COMM_WORLD);


  ri = rank * my_nrows; //-1
  rf = ri + my_nrows;//+1
  posi = ri*n;
  posf = rf*n;

 while ( error > tol*tol && iter < iter_max )
  {
    iter++;
    //MPI_Send(buffer, count , type ,dest, tag, comm);
    if(rank > MASTER){
      MPI_Send(my_A+n, n, MPI_FLOAT, rank-1, tag ,MPI_COMM_WORLD);
      MPI_Recv(my_A  , n, MPI_FLOAT, rank-1, tag ,MPI_COMM_WORLD, &Stat);
    }
    if(rank < numtasks -1 ){
       MPI_Send(  (my_A + n*(my_nrows))  , n, MPI_FLOAT, rank+1, tag ,MPI_COMM_WORLD);
       MPI_Recv( (my_A + n*(my_nrows+1) )  , n, MPI_FLOAT, rank+1, tag ,MPI_COMM_WORLD, &Stat);
    }

    /*MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 1){
    printf("La matriu del procés %d és: \n", rank);
    print_matrix(my_A, my_nrows+2,n);
    printf("\n");
    }
    

    MPI_Barrier(MPI_COMM_WORLD);

    exit(0);
    if(rank == MASTER){
        printf("He passat la barrera (Master)");
    }
    */


    my_error= my_laplace_step(my_A, my_temp, my_nrows +2, n);

    MPI_Reduce(&my_error, &error, 1, MPI_FLOAT, MPI_MAX, MASTER, MPI_COMM_WORLD);

    //Send the new portions of matrix to MASTER
    /*MPI_Gather(
    void* send_data,
    int send_count,
    MPI_Datatype send_datatype,
    void* recv_data,
    int recv_count,
    MPI_Datatype recv_datatype,
    int root,
    MPI_Comm communicator)*/





    /*if(rank != MASTER){
          MPI_Send(&my_A[n] , n*my_nrows, MPI_FLOAT, MASTER, tag ,MPI_COMM_WORLD);
    }

    if(rank == MASTER){
        for(int id = 1; id< numtasks;id++){
          //assuming my_nrows of the MASTER is the same as for the others
            MPI_Recv(&A[id*my_nrows*n] , n*my_nrows, MPI_FLOAT, id, tag ,MPI_COMM_WORLD, &Stat);
       }
       my_init(my_A,A, 1 , my_nrows-2,ri, rf, n); 
    }
    */



    //float *swap= A; A=temp; temp= swap; // swap pointers A & temp
    float *swap= my_A; my_A=my_temp; my_temp= swap;

    //Bcast of error again?
  }





 //}


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
