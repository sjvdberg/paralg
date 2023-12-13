#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    int p, s;

    /* SPMD part */
    MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD,&p); /* p = number of processors */ 
    MPI_Comm_rank(MPI_COMM_WORLD,&s); /* s = processor number */ 

    int n=1000;
    if (argc>1) n=atoi(argv[1]);

    if(n<0) MPI_Abort(MPI_COMM_WORLD,-1);

    computeVector(n, p, s, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
} /* end main */

//n is number of rows/columns.
//p is number of processors.
//s is own rank.
void computeVector(int N, int p, int s, MPI_Comm comm)
{

    int numrows = (N+p-s-1)/p ;
    int remainder = N % p;
    int firstrow;
    if( s < remainder) 
        //numrows is the same as all earlier processors.
        firstrow = s * numrows;
    else
        //The first remainder processors have 1 more element
        firstrow = s * numrows + remainder;
    
    int baseRows[numrows][11];
    for(int i = 0; i < numrows; i++)
    {
        int k = rand() % 10;
        for(int l = 0; l < 11; l++)
        {
            if(l <= k)
                baseRows[i][l] = rand() % N;
            else
                baseRows[i][l] = -1;
        }
    }
    printf("proces %i has computed their rows\n", s);
    int localDiagonal[N];
    for(int i = 0; i < N; i++)
        localDiagonal[i] = 0;
    for(int i = 0; i < numrows; i++)
        for(int l = 0; l < 10; l++)
            if(baseRows[i][l] != -1)
                localDiagonal[baseRows[i][l]]++;
    int Diagonal[N];
    for(int i = 0; i < N; i++)
        printf("processor %i. %i\n", s, localDiagonal[i]);
    MPI_Request requests[2*p];
    for(int r = 0; r < p; r++)
    {
        if(r == s)
        {
            for(int i = 0; i < numrows; i++)
                Diagonal[i] = localDiagonal[i];
            continue;
        }
        MPI_Isend(localDiagonal, 2*N, MPI_INT, r, r, comm, &requests[r]);
    }
    MPI_Barrier(comm);
    for(int r = 0; r < p; r++)
    {
        if(r == s) continue;
        MPI_Irecv(localDiagonal, 2*N, MPI_INT, r, s, comm, &requests[p+r]);
        for(int i = 0; i < numrows; i++)
            Diagonal[i] += localDiagonal[i];
    }
    printf("succesfully computed diagonal.\n");
    for(int i = 0; i < N; i++)
        printf("processor %i. %i . %i\n", s, localDiagonal[i], Diagonal[i]);
}