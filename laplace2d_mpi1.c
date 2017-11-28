/*  Copyrigh (c) 2012, NVIDIA CORPORATION. All rights reserved.
 *
 * redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 *  are met:
 *   * Redistributions of source code must retain the above copyright
 *   *    notice, this list of conditions and the following disclaimer.
 *   *  * Redistributions in binary form must reproduce the above copyright
 *   *    notice, this list of conditions and the following disclaimer in the
 *   *    documentation and/or other materials provided with the distribution.
 *   *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *   *    contributors may be used to endorse or promote products derived
 *   *    from this software without specific prior written permission.
 *   *
 *   * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 *   * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *   * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *   * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 *   * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *   * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *   * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *   * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 *   * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *   */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#ifndef n
#define n 4096
#endif
#ifndef m
#define m 4096
#endif

float A[n][m];
float Anew[n][m];
float y[n];

int main(int argc, char** argv)
{
    int i, j;
    int iter_max = 1000;
    
    const float pi  = 2.0f * asinf(1.0f);
    const float tol = 3.0e-3f;
    float error= 1.0f;    

    // obtener argumentos proporcionados en tiempo de ejecucion
    if (argc>1) {  iter_max = atoi(argv[1]); }

    memset(A, 0, n * m * sizeof(float));
    
    //  set boundary conditions
    for (i=0; i < m; i++)
    {
       A[0][i]   = 0.f;
       A[n-1][i] = 0.f;
    }

    for (j=0; j < n; j++)
    {
       y[j] = sinf(pi * j / (n-1));
       A[j][0] = y[j];
       A[j][m-1] = y[j]*expf(-pi);
    }

    printf("Jacobi relaxation Calculation: %d x %d mesh, maximum of %d iterations\n", 
           n, m, iter_max );

    int iter = 0;

    int numtasks, rank, dest, source, tag = 1,rc;

    rc = MPI_Init (&argc, &argv);
    if (rc != MPI_SUCCESS)
    {
      printf ("Error starting MPI program. Terminating.\n");
      MPI_Abort (MPI_COMM_WORLD, rc);
      return -1;
    }

    MPI_Comm_size (MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    if(rank == 0){

    while ( error > tol && iter < iter_max )
    {
       for( i=1; i < m-1; i++ )
          for( j=1; j < n-1; j++)
              Anew[j][i] = ( A[j][i+1]+A[j][i-1]+A[j-1][i]+A[j+1][i]) / 4;

       error = 0.f;
       for( i=1; i < m-1; i++ )
          for( j=1; j < n-1; j++)
              error = fmaxf( error, sqrtf( fabsf( Anew[j][i]-A[j][i] ) ) );

       for( i=1; i < m-1; i++ )
          for( j=1; j < n-1; j++)
               A[j][i] = Anew[j][i];

       iter++;
       if(iter % (iter_max/10) == 0) printf("%5d, %0.6f\n", iter, error);
    }

  }
    MPI_Finalize();

    return 0;
}
