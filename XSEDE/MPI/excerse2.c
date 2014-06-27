#include <stdio.h>
#include <mpi.h>

int main( int argc, char **argv ){

 int max_pe_num, my_pe_num, numprocs, index;

 MPI_Init(&argc, &argv);
 MPI_Comm_rank(MPI_COMM_WORLD, &my_pe_num );

 MPI_Reduce(&my_pe_num, &max_pe_num, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);

 if( my_pe_num == 0 ){
 printf("Number of procs: %d.\n", max_pe_num);
 }

 MPI_Finalize();
}

