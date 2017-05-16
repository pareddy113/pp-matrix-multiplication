/*
 Approach2: Cannon Algorithm
 Cannons algorithm is a distributed algorithm for Matrix multiplication for a two dimensional meshes of
 size sqrt(p) × sqrt(p) processors. The storage space remains same irrespective of number of processors
 Topology:
 2-D Mesh Topology with wrap around are used to implement this algorithm. Ring is used for smaller P.
 The design of the algorithm is discussed in further sections. The reason for using a mesh topology is to
 exploit the number of connected nodes to ease communication. This algorithm is more suited for mesh
 than any other topology due to the circular left and up shifts after each iteration.
 Method 1:
 Program Cannon_finalMPI.c
 Job finalcannon.pbs
 Output file finalcannon.out
 In this method, root processor creates the array of size N * N and populates it with random values. Input
 method one is chosen for populating with values in range {-1, 0, 1}. When N = sqrt(p), each processor has
 elements aij, bij which correspond to the matrices A, B which are to be multiplied. If N > sqrt(p), blocks of
 rows N/sqrt(p)
 could be assigned to each processor.
 The main idea behind this algorithm is to perform sqrt(p) × sqrt(p) multiplications and additions in parallel and
 then supplement the data needed by talking to neighbors for data (4N-2 messages);
 There are three main calculations to be made when this algorithm is used.
 1. Skew step :
 a. Each row of A[i] [] is circularly left shifted by j units wherei+ j is the column number
 b. Each column of B[][j] is circularly shifted up by i+j units where I is the row number
 2. Computation
 a. Each element A[i][j] and B[i][j] are multiplied at the processor Pi,j and the partial sum is
 stored in C[i][j]
 3. Circular shift
 After computing the partial C[i][j], next data at Pi,j needs to be fetched according the the
 following rule
 a. Each row of A[i][j] is circularly left shifted by 1 unit (take data from right neighbor and
 pass current data to left one)
 b. Each column of B[i][j] is circularly up shifted by 1 unit. (Take data from bottom neighbor
 and pass current to top neighbor)
 4. Above three steps result in calculating A X A, to calculate A^(k), the above 3 steps are repeated by
 making A = C (result from A X A) and retaining B = A. This results in the above code being in a
 loop for k-1 iterations
 5. After A^(k) is calculated, LU decomposition is used to calculate the determinant of A^(k).
 a. Division step: Each element of A[i][k] of kth row is divided by A[k][k]
 b. Communication: Each element in the kth row is needed in further rows for elimination
 hence these are broadcasted to the rows below the kth row
 c. Elimination: Elimination of uses the following equation A[I,j] = A[I,j] – A[I,k]*A[i][k]
 d. Gather all the diagonal eliments of U inLU and multiply them to get the determinant
 For the communication it can be seen that the following strategy of 2D-Mesh would make the
 communication efficient
 */


/*
This program does A^k of a matrix a using cannons algorithm on p processors
*/
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
/***************************************************************************************
//functions 
****************************************************************************************/
/***************************************************************************************
//Function creates  a matrix of given size 
****************************************************************************************/
float** intialize(int size)
{
	int i, j;
	i = j = 0;
	float *data = (float *)malloc(size*size*sizeof(float));
	float** rand_nums = (float **)malloc(sizeof(float *)* size);
	for (i = 0; i < size; i++) {
		rand_nums[i] = &(data[size*i]) ;
		//rand_nums[i] = (float *)malloc(sizeof(float)*size);
	}
	return rand_nums;	
}
/***************************************************************************************
//Function populates a matrix with 1,0-1 radomly based on probability of generating 0,1,2
****************************************************************************************/
float** randomizeMatrix(int size)
{
	int i, j,count;
	i = j = count=0;
	float *data = (float *)malloc(size*size*sizeof(float));
	float** rand_nums = (float **)malloc(sizeof(float *)* size);
	for (i = 0; i < size; i++) {
		rand_nums[i] = &(data[size*i]) ;
	}
  	for (i = 0; i < size; i++) {
		for (j =0; j<size;j++)
		{
			int k = rand()%3;
			int z = 0;
			count++;
			switch(k)
			{
				case 1: z = -1;
					break;
				case 2: z = 1;
					break;
				default: z=0;
			}
			rand_nums[i][j] = z;
//			rand_nums[i][j] = count;
			if(size<32)
			printf("%f ",rand_nums[i][j]); 
  		}
		if(size<32)
  		printf("\n");
	}
	return rand_nums;	
}
/***************************************************************************************
//Function does a circular upward shift by 1 of matrix 
****************************************************************************************/
float** circular_up(float** matrix, int size,int per_proc, MPI_Comm ring, int rank, int procs)//, int i, int j)
{
	int i,j;
	float** temp = intialize(size);
	int recvfrom = (rank+1)%procs;
	int sendto = (rank-1)%procs;
	MPI_Status status;
	if(sendto<0)
	{
		sendto = procs-1;
	}
//	printf("sendto %d, rec %d\n procs %d",sendto,recvfrom,procs);
	for (i = 0; i < per_proc-1; i++) {
		for (j =0; j<size;j++)
		{
		//	int k = (i + 1)%size;
			temp[i][j] = matrix[i+1][j];
		}
	}
	float* rotate = (float *)malloc(size*sizeof(float));
	if(rank%2==0)
	{
		MPI_Send(&(matrix[0][0]),size,MPI_FLOAT,sendto,rank,ring);
		MPI_Recv(rotate,size,MPI_FLOAT,recvfrom,recvfrom,ring,&status);
	}
	else
	{
		MPI_Recv(rotate,size,MPI_FLOAT,recvfrom,recvfrom,ring,&status);
		MPI_Send(&(matrix[0][0]),size,MPI_FLOAT,sendto,rank,ring);
	}
	for(j =0;j<size;j++)
	{
		temp[per_proc-1][j] = rotate[j];
	}
//	return temp;
/*	if (rank == 0)
	{
	printf("\ncircular up %d\n",rank);
	for (i =0;i<per_proc;i++)
	{
		for(j=0;j<size;j++)
		{
			printf("%f ",temp[i][j]);
		}
		printf("\n");
	}
	}*/
	return temp;
}
/***************************************************************************************
//Function does a circular left shift by 1 of matrix 
****************************************************************************************/
float** circular_left(float** matrix, int size,int per_proc,int rank)//, int i, int j)
{
	int i,j;
	float** temp = intialize(size);
	for (i = 0; i < per_proc; i++) {
		for (j =0; j<size;j++)
		{
			int k = (j + 1)%size;
			temp[i][j] = matrix[i][k];
		}
	}
/*	if(rank == 0)
	{
	printf("\ncircular left \n");
	for (i =0;i<per_proc;i++)
	{
		for(j=0;j<size;j++)
		{
			printf("%f ",temp[i][j]);
		}
		printf("\n");
	}
	}*/
	return temp;
}
/***************************************************************************************
//Function does a skew left of original matrix before sending
****************************************************************************************/
float** skew_left(float** matrix,int size)
{
	int i,j;
	float** temp = intialize(size);
	for (i = 0; i < size; i++) {
		for (j =0; j<size;j++)
		{
			int k = (j + i)%size;
			temp[i][j] = matrix[i][k];
		}
	}
	/*printf("\n skew left\n");

	for (i =0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			printf("%f ",temp[i][j]);
		}
		printf("\n");
	}
	//	printf("address of tem = %d\n",temp);
*/
	return temp;
}
/***************************************************************************************
//Function does a skew up of original matrix before sending
****************************************************************************************/
float** skew_up(float** matrix,int size)
{
	int i,j;
	float** temp = intialize(size);
	//printf("address of tem = %d\n",temp);

	for (i = 0; i < size; i++) {
		for (j =0; j<size;j++)
		{
			int k = (i + j)%size;	
	//		printf("i = %d, j=%d, k=%d ",i,j,k);
			temp[i][j] = matrix[k][j];
		//		printf(" j = %d",j);
	//		printf("matrix[%d][%d] = %f, temp[%d][%d]=%f \n",k,j,matrix[k][j],i,j,temp[i][j]);
		}
	}
	//printf("address of tem = %d\n",temp);

/*	printf("\nskew up\n");
	for (i =0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			printf("%f ",temp[i][j]);
		}
		printf("\n");
	}*/
	return temp;
}
/***************************************************************************************
//Function does the main matrix multiplication using cannons method
****************************************************************************************/
float** cannon(float** matrix,float** matrixB,int size,int procs,int per_proc,int rank,MPI_Comm ring)//, MPI_Comm comm)
{
	int i,j,k;
	float** temp = intialize(size);
//	printf("\nproduct %d\n", rank);
	for(i =0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			temp[i][j] = 0.0;
		}
	}
//	printf(" start main loop\n");
	for(k=0;k<size;k++)
	{
		for (i =0;i<per_proc;i++)
		{
			for(j=0;j<size;j++)
			{
				temp[i][j] = temp[i][j] + matrix[i][j] * matrixB[i][j];
			}
		}
//		printf("before circular stuff\n");
		matrix = circular_left(matrix,size,per_proc,rank);
//		printf("before matb circular up%d\n",procs);
		matrixB = circular_up(matrixB,size,per_proc,ring,rank,procs);
//		printf("after circular stuff\n"); 
/*		int count =0;
		if(rank == 7 && k == 0)
	        {
			printf("inside print \n");
			for(i =0;i<per_proc;i++)
			{
				for(j=0;j<size;j++)
				{	
					count++;
					printf("%f ",matrixB[i][j]);
				}
				printf("\n");
			}
		printf("/n circular up %d rank %d\n",count,rank);
		}*/
	//	printf("/n circular up %d rank %d\n",count,rank);
//float** circular_up(float** matrix, int size,int per_proc, MPI_Comm ring, int rank, int procs)//, int i, int j)
	}
  //     printf(" before gather %d rank %d per\n",rank,per_proc);
	// gather all products at root
       //MPI_Gather(&(temp[0][0]),size*per_proc,MPI_FLOAT,&(temp[rank*per_proc][0]),size*per_proc,MPI_FLOAT,0,ring);
	int recv_size = per_proc*size;
	MPI_Status  status;
	if(rank == 0)
	{
		for(i =1;i<procs;i++)
		{
			MPI_Recv(&(temp[i*per_proc][0]),recv_size,MPI_FLOAT,i,1,ring,&status);
		}
	}
	else
	{
		MPI_Send(&(temp[0][0]),recv_size,MPI_FLOAT,0,1,ring);
	}
//	printf("after gather%d rank %d per\n",rank,per_proc);
	return temp;

}
//send matrix from root;
float** send_matrix(float** matrix, int size, MPI_Comm ring, int procs, int per_proc,int rank);
float find_det(float** matrix, int size);//,int procs, int per_proc,int rank, MPI_Comm ring);
/* Main processes, processing starts here*/
int main (int argc, char** argv)
{
	int power = 5;
	int size = 4;
//	int procs = 2;
	//int per_proc = size/procs;
	float** matrix = intialize(size);
	float** matrixB = intialize(size);
	float** multMatrix = intialize(size);
	//Declare the MPI variables
	int rank,iter;
 	int world_rank,world_size;
	int dim =1;
 	int periods =1;
 	double start,end,send,decomp,exponent,comm,setup,accum;
	accum = 0.0;
 	MPI_Init(&argc,&argv);
 	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
 	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 	int procs = world_size;
	int per_proc = size/procs;
//	printf("%d world size %d per proc %d size \n",procs,per_proc,size);
 	//create a ring of size processors
	int process_per_dim =world_size ;
	MPI_Comm ring;
	MPI_Cart_create(MPI_COMM_WORLD, dim,&process_per_dim,&periods,0,&ring);
	MPI_Comm_rank(ring,&rank);
	MPI_Status status;
	MPI_Barrier(ring);
	//populate array at rank 0
	start =MPI_Wtime();
	if(rank==0)
	{
		float** temp = intialize(size);
		temp = randomizeMatrix(size);
		//intial skew of the matrix and send it to all procs
		matrix = skew_left(temp,size);
		matrixB = skew_up(temp,size);
		free (temp);
	}
	MPI_Barrier(ring);
	end=MPI_Wtime();
	setup = end-start;
	accum +=setup;
//	printf(" setuptiime = %f \n",setup);
//	printf("begining exponentiation\n");
	start = MPI_Wtime();
	for(iter =1; iter<power;iter++)
	{
                double start1 = MPI_Wtime();		
		matrix = send_matrix(matrix,size,ring,procs,per_proc,rank);
		matrixB = send_matrix(matrixB,size,ring,procs,per_proc,rank);
		MPI_Barrier(ring);
		double end1 = MPI_Wtime();
		send = end1-start1;
	/*printf("\ncircular left\n");
	multMatrix = circular_left(matrix,size);
	printf("\ncricular up\n");
	multMatrix = circular_up(matrix,size);*/
	//MPI_Barrier(ring);
//	printf("per_proc %d rank %d \n",per_proc,rank);
		int i,j,count;
		count =0;
//		MPI_Barrier(ring);
/*		if(rank == 3)
		{
	//printf("inside print \n");
			for(i =0;i<per_proc;i++)
			{
				for(j=0;j<size;j++)
				{	
					count++;
	//		printf("%f ",matrixB[i][j]);
			}
	//	printf("\n");
			}
		}		*/
//	printf("/n cannon product %d rank %d\n",count,rank);
		start1 = MPI_Wtime();
		multMatrix = cannon(matrix, matrixB,size,procs,per_proc,rank,ring);
		MPI_Barrier(ring);
		end1=MPI_Wtime();
		comm = end1-start1;
//	printf("finished product %d\n",rank);
		start1 = MPI_Wtime();
		if(rank == 0)
		{
			float** temp = intialize(size);
			printf("inside print prod\n");
			for(i =0;i<size;i++)
			{	
				for(j=0;j<size;j++)
				{	
	//			count++;
				temp[i][j] = multMatrix[i][j];
		//		printf("%f ",multMatrix[i][j]);
				}
		//	printf("\n");
			}
			if(iter+1<power)
			{
			matrix = skew_left(temp,size);
			}
		}
		MPI_Barrier(ring);
		end1 = MPI_Wtime();
		decomp = end1-start1;
		if(rank==0)
		{      
			accum += decomp+comm+send;
			printf("\n power %d time %f procs %d size %d\n",iter+1,accum, procs,size);
		}
	}
	end = MPI_Wtime();
	exponent = end-start;
	//send multMatrix to calculate det
	//matrix = send_matrix(multMatrix,size,ring,procs,per_proc,rank);
	float det = 0.0;
	if(rank == 0)
	{det = find_det(multMatrix,size);//,procs,per_proc,rank,ring);
	 printf("det = %f, %f\n",det,(exponent+setup) );
	}
//float deta(float** matrix, int size,int procs, int per_proc,int rank, MPI_Comm ring);
//	printf("/n cannon product %d rank %d\n",count,rank);
//	free(matrix);
//	free(matrixB);
//	free(multMatrix);//*/
	MPI_Finalize();
	return 0;
}
//send matrix across 
float** send_matrix(float** matrix, int size, MPI_Comm ring, int procs, int per_proc,int rank)
{
	int recv_size = per_proc*size;
	MPI_Status status;
	int i;
//	printf("starting send%d %d\n", rank,procs);
	if(rank == 0)
	{
		for(i=1;i<procs;i++)
		{
			MPI_Send(&(matrix[i*per_proc][0]),recv_size,MPI_FLOAT,i,1,ring);
		}
//	printf("starting recv%d\n",rank);
	}
	//printf("starting recv\n");
	else//(rank > 0)
	{
		MPI_Recv(&(matrix[0][0]),recv_size,MPI_FLOAT,0,1,ring,&status);//MPI_STATUS_IGNORE);//&status);
//		printf("recv over for %d, first elm is %f \n",rank,matrix[0][0]);
	}
	return matrix;
} 
float find_det (float **array, int size)
{ 
	int i,j,k;
	for (k=0;k<size;k++)
	{
	//	printf("inside k %d\n",k);
		for(j=k+1;j<size;j++)
		{
			array[j][k]= (array[j][k])/(array[k][k]);

		}
	//	printf("after first j");
		for (j=k+1 ; j<size;j++)
		{
			for(i= k+1;i<size;i++)
			{
				array[i][j] = array[i][j]- (array[i][k]*array[k][j]);
			}
		}
	} 
	//printf("multiply \n");
	//determinant using U
	float det =1;
	for (k=0;k<size;k++)
	{
	/*	for (j=k;j<size;j++)
		{
                   printf(" %f ",array[k][j]);
		}*/
		det = det * array[k][k];
	//	printf("\n");
	}//*/
        return det;
	
}


		




