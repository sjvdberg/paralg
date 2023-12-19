#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int process_Rank, size_Of_Cluster, message_Item;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

    if(process_Rank == 0){
        int arr[20];
        for(int i = 0; i < 20; i++)
            arr[i] = i;
        
        MPI_Send(&arr[0], 5, MPI_INT, 1, 1, MPI_COMM_WORLD);
        printf("Message Sent:\n");
        for(int i = 0; i < 5; i++)
            printf("%i. %i\n", process_Rank, i);
    }

    else if(process_Rank == 1){
        int receive[5];
        MPI_Recv(receive, 5, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Message Received:\n");
        for(int i = 0; i < 5; i++)
            printf("%i. %i\n", process_Rank, receive[i]);
    }

    MPI_Finalize();
    return 0;
}