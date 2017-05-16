/*
 1D Algorithm
 
 In this algorithm,
 1. the rows of the matrices are divided between processes, so each row is sent to each process.
 2. Each process receives its share of matrices and performs multiplication with all the columns to
 get the A square.
 3. Calculated output is sent back to the source node.
 4. Then based on the k value, the source node sends it back to the worker processes until k value is
 reached.
 5. The final output is sent back to the source node.
 6. Then the determinant is calculated using the serial determinant algorithm.
 
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define matrix_size 128                 
#define origin_master 1          
#define origin_worker 2          

float find_det (float **array, int size)
{ 
 int i,j,k;
 for (k=0;k<size;k++)
 {
  for(j=k+1;j<size;j++)
  {
   array[j][k]= (array[j][k])/(array[k][k]);
  }
  for (j=k+1 ; j<size;j++)
  {
   for(i= k+1;i<size;i++)
   {
    array[i][j] = array[i][j]- (array[i][k]*array[k][j]);
   }
  }
 } 
 float det =1;
 for (k=0;k<size;k++)
 {
  for (j=k;j<size;j++)
  {
                   printf(" %f ",array[k][j]);
  }
  det = det * array[k][k];
  printf("\n");
 }
        return det;
 
}

int main (int argc, char *argv[]) {
int world_size,processID,numworkers,source,dest,mes_type,rows,rows_per_process, extra, offset,i, j,k,rc,t=0,m=0,n=0,MASTER=0,power=6,v;
double total_time,total_comm_time,total_comp_time;
double a[matrix_size][matrix_size],b[matrix_size][matrix_size],c[matrix_size][matrix_size],d[matrix_size][matrix_size];

//int **c=(int **)malloc(matrix_size*sizeof(int));
//    for(i=0;i<matrix_size;i++)
  //      c[i]=(int *)malloc(matrix_size*sizeof(int));

total_time=total_comm_time=total_comp_time=0.0;

MPI_Status status;
MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&processID);
MPI_Comm_size(MPI_COMM_WORLD,&world_size);
if (world_size < 2 ) {
  printf("Need at least two MPI tass. Quitting...\n");
  MPI_Abort(MPI_COMM_WORLD, rc);
  exit(1);
  }
numworkers = world_size-1;

total_time=total_time-MPI_Wtime();
//master process

   total_comp_time-=MPI_Wtime();
   if (processID == 0)
  {
   printf("Total number of processors are %d \n",world_size);
      printf("Array A is =\n");
	
      for (i=0; i<matrix_size; i++)
      {
         for (j=0; j<matrix_size; j++)
       {
            v=rand();
	switch(v%3)
	{
	    case 1:
		a[i][j]=1;
		break;
	    case 2:
		a[i][j]=-1;
		break;
	    default:
		a[i][j]=0;
         }
        }
       }

	for (i=0; i<matrix_size; i++)
	{
          for (j=0; j<matrix_size; j++)
 	 {
           printf("%6.2f ", a[i][j]);
	 }
	printf("\n");
	}

	total_comp_time+=MPI_Wtime();
      /* Divide the work among worker process*/
      rows_per_process = matrix_size/numworkers;
      extra = matrix_size%numworkers;
      offset = 0;
      mes_type = origin_master;
       
      for (dest=1; dest<=numworkers; dest++)
      {
	if (dest<=extra)
		rows=rows_per_process+1;
	else
		rows=rows_per_process;
	
	total_comm_time-=MPI_Wtime();
         MPI_Send(&offset, 1, MPI_INT, dest, mes_type, MPI_COMM_WORLD);
         MPI_Send(&power, 1, MPI_INT, dest, mes_type, MPI_COMM_WORLD);
         MPI_Send(&rows, 1, MPI_INT, dest, mes_type, MPI_COMM_WORLD);
         MPI_Send(&a[offset][0], rows*matrix_size, MPI_DOUBLE, dest, mes_type,
                   MPI_COMM_WORLD);
         MPI_Send(&a, matrix_size*matrix_size, MPI_DOUBLE, dest, mes_type, MPI_COMM_WORLD);
	
	total_comm_time+=MPI_Wtime();
         offset = offset + rows;
      }
     
      /* Receive computed results from worker process*/
      mes_type = origin_worker;
      for (i=1; i<=numworkers; i++)
      {
         source = i;

	total_comm_time-=MPI_Wtime();
         MPI_Recv(&offset, 1, MPI_INT, source, mes_type, MPI_COMM_WORLD, &status);
         MPI_Recv(&rows, 1, MPI_INT, source, mes_type, MPI_COMM_WORLD, &status);
         MPI_Recv(&c[offset][0], rows*matrix_size, MPI_DOUBLE, source, mes_type, 
                  MPI_COMM_WORLD, &status);
	
	total_comm_time+=MPI_Wtime();
      }
	
	total_comp_time-=MPI_Wtime();
      /* Print the final matrix */
      printf("******************************************************\n");
      printf("Result Matrix:\n");
      for (i=0; i<matrix_size; i++)
      {
         printf("\n"); 
         for (j=0; j<matrix_size; j++) 
            printf("%6.2f   ", c[i][j]);
      }
      printf("\n******************************************************\n");
      printf ("Done.\n");
//	printf("The determinant is= %f",find_det(c,matrix_size));	
	total_comp_time+=MPI_Wtime();
   }

//worker process
   if(processID > 0)
   {
      mes_type = origin_master;
	
	total_comm_time-=MPI_Wtime();
      MPI_Recv(&offset, 1, MPI_INT, MASTER, mes_type, MPI_COMM_WORLD, &status);
      MPI_Recv(&power, 1, MPI_INT, MASTER, mes_type, MPI_COMM_WORLD, &status);
      MPI_Recv(&rows, 1, MPI_INT, MASTER, mes_type, MPI_COMM_WORLD, &status);
      MPI_Recv(&a, rows*matrix_size, MPI_DOUBLE, MASTER, mes_type, MPI_COMM_WORLD, &status);
      MPI_Recv(&b, matrix_size*matrix_size, MPI_DOUBLE, MASTER, mes_type, MPI_COMM_WORLD, &status);
	
	total_comm_time+=MPI_Wtime();

	total_comp_time-=MPI_Wtime();
    for(t=0;t<power;t++)
     {
       for(m=0; m<rows;m++){
	for(n=0;n<matrix_size;n++){
	 d[m][n]= a[m][n];
//	printf("%d--%d--%d--%6.2f \n",processID,m,n,d[m][n]);
       }
      }

      for (k=0; k<matrix_size; k++)
       {
         for (i=0; i<rows; i++)
         {
            a[i][k] = 0.0;
            for (j=0; j<matrix_size; j++)
	     {
               a[i][k] = a[i][k] + d[i][j] * b[j][k];
	     }
         }
       }
     }
	
	total_comp_time+=MPI_Wtime();
	total_time+=MPI_Wtime();	
//	total_comp_time+=MPI_Wtime();
      mes_type = origin_worker;
      MPI_Send(&offset, 1, MPI_INT, MASTER, mes_type, MPI_COMM_WORLD);
      MPI_Send(&rows, 1, MPI_INT, MASTER, mes_type, MPI_COMM_WORLD);
      MPI_Send(&d, rows*matrix_size, MPI_DOUBLE, MASTER, mes_type, MPI_COMM_WORLD);
}
	MPI_Barrier(MPI_COMM_WORLD);
    if(processID==0){

    printf("\nThe total communication time = %f seconds\n",total_comm_time);
    printf("\nThe total computation time = %f seconds\n",total_comp_time);
    printf("\nThe total running time = %f seconds\n",total_comm_time+total_comp_time);

	}
    MPI_Finalize();
}
