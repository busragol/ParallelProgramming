/////////////////////////////////////////////////////////////////////
// CENG316 Parallel Programming
// Programing exam 1
//
// Author: Group???
// Group members:  Aaa Bbbb , Cccc Dddd
// Date: 24/03/2019
// Description: implementation of myMPI_Allreduce
/////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h> 

#define NUMBER_OF_TESTS 1000
#define ARRAY_SIZE 100000

void myMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int main(int argc, char **argv )
{
    int          rank,  numproc;
    double       t1, t2;
    int		     d_in[ARRAY_SIZE];
    int 		 d_out_sum[ARRAY_SIZE], d_out_max[ARRAY_SIZE],d_out_min[ARRAY_SIZE];
    int          i, j, k, nloop;
	
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &numproc );
    srand(rank+1); 
    
    if (rank == 0) 
		printf( "Test launched with %d processor, ARRAY_SIZE: %d \n",numproc, ARRAY_SIZE );
    
    for(i=0;i<ARRAY_SIZE;i++)
		d_in[i] = rand()% 100;
		
	MPI_Barrier(MPI_COMM_WORLD );			
	t1 = MPI_Wtime();		
    for (k=0; k<NUMBER_OF_TESTS; k++) {		
		MPI_Allreduce( d_in, d_out_sum, ARRAY_SIZE, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
		MPI_Allreduce( d_in, d_out_max, ARRAY_SIZE, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
		MPI_Allreduce( d_in, d_out_min, ARRAY_SIZE, MPI_INT, MPI_MIN, MPI_COMM_WORLD );		
    }
    t2 = MPI_Wtime() - t1;		
    	
    MPI_Barrier(MPI_COMM_WORLD);
    int avg=0;
	for(i=0;i<ARRAY_SIZE;i++){ avg += d_out_sum[i]; } avg = avg / (numproc*ARRAY_SIZE);
	printf("Proc %2d: Average of all elements: %d \n", rank, avg );
	
    if (rank == 0) {		
		printf( "\nAverage time for Allreduce: %.5f msec\n", t2/NUMBER_OF_TESTS );
		printf( "First 3-elements of sum array: %d %d %d \n", d_out_sum[0], d_out_sum[1],d_out_sum[2]);		
		printf( "First 3-elements of max array: %d %d %d \n", d_out_max[0], d_out_max[1], d_out_max[2] );		
		printf( "First 3-elements of min array: %d %d %d \n", d_out_min[0], d_out_min[1], d_out_min[2] );		
		printf( "Range of First 3-elements: \t%d %d %d \n", d_out_max[0]-d_out_min[0], d_out_max[1]-d_out_min[1], d_out_max[2]-d_out_min[2]  );
		printf( "Average of First 3-elements: \t%d %d %d \n\n",  d_out_sum[0]/numproc, d_out_sum[1]/numproc,d_out_sum[2]/numproc );		
    }
    
    // reinitializing output arrays
    for(i=0;i<ARRAY_SIZE;i++){ d_out_sum[i] = 0; d_out_max[i]=0; d_out_min[i]=0;  } 
    
    ///////////////////////////////////////////////////////////////////
    /* Calling your implementation of AllReduce */
    /* myMPI_Allreduce */
    
    MPI_Barrier(MPI_COMM_WORLD );			
	t1 = MPI_Wtime();
    for (k=0; k<NUMBER_OF_TESTS; k++) {				
		myMPI_Allreduce( d_in, d_out_sum, ARRAY_SIZE, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
		myMPI_Allreduce( d_in, d_out_max, ARRAY_SIZE, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
		myMPI_Allreduce( d_in, d_out_min, ARRAY_SIZE, MPI_INT, MPI_MIN, MPI_COMM_WORLD );		
    }
    t2 = MPI_Wtime() - t1;	
    
	MPI_Barrier(MPI_COMM_WORLD);
    avg=0;
	for(i=0;i<ARRAY_SIZE;i++){ avg += d_out_sum[i]; } avg = avg / (numproc*ARRAY_SIZE);
	printf("Proc %2d: Average of all elements: %d \n", rank, avg );	
    		
    if (rank == 0) {
		printf( "\nAverage time for myAllreduce: %.5f msec\n", t2/NUMBER_OF_TESTS );		
		printf( "First 3-elements of sum array: %d %d %d \n", d_out_sum[0], d_out_sum[1],d_out_sum[2]);		
		printf( "First 3-elements of max array: %d %d %d \n", d_out_max[0], d_out_max[1], d_out_max[2] );		
		printf( "First 3-elements of min array: %d %d %d \n", d_out_min[0], d_out_min[1], d_out_min[2] );		
		printf( "Range of First 3-elements: \t%d %d %d \n", d_out_max[0]-d_out_min[0], d_out_max[1]-d_out_min[1], d_out_max[2]-d_out_min[2]  );
		printf( "Average of First 3-elements: \t%d %d %d \n\n",  d_out_sum[0]/numproc, d_out_sum[1]/numproc,d_out_sum[2]/numproc );
    }

    MPI_Finalize( );
    return 0;
}

void myMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm){
// write your code here
	
}
