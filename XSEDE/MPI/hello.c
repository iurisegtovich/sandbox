#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv){
    int my_PE_num;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_PE_num);
    printf("Hello World, PE num %d.\n",my_PE_num);
    MPI_Finalize();
    return 0;
}
