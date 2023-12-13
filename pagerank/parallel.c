#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

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


int numRows(int N, int p, int s)
{
    return (N+p-s-1)/p;
}
int firstRow(int N, int p, int s)
{
    int numrows = numRows(N, p,s);
    int remainder = N % p;
    int firstrow;
    if( s < remainder) 
        //numrows is the same as all earlier processors.
        firstrow = s * numrows;
    else
        //The first remainder processors have 1 more element
        firstrow = s * numrows + remainder;
    return firstrow;
}

//n is number of rows/columns.
//p is number of processors.
//s is own rank.
void computeVector(int N, int p, int s, MPI_Comm comm)
{
    srand(time(0));
    int numrows = (N+p-s-1)/p ;
    int firstrow = firstRow(N, p, s);
    
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
    int numElements = 0;
    for(int i = 0; i < N; i++)
        localDiagonal[i] = 0;
    for(int i = 0; i < numrows; i++)
        for(int l = 0; l < 10; l++)
            if(baseRows[i][l] != -1)
            {
                numElements++;
                localDiagonal[baseRows[i][l]]++;
            }
    float Diagonal[N];
    MPI_Request requests[2*p];
    for(int r = 0; r < p; r++)
    {
        if(r == s)
        {
            for(int i = 0; i < N; i++)
                Diagonal[i] = localDiagonal[i];
            continue;
        }
        MPI_Isend(localDiagonal, N, MPI_INT, r, r, comm, &requests[r]);
    }
    MPI_Barrier(comm);
    for(int r = 0; r < p; r++)
    {
        if(r == s) continue;
        MPI_Irecv(localDiagonal, N, MPI_INT, r, s, comm, &requests[p+r]);
        for(int i = 0; i < N; i++)
            Diagonal[i] += localDiagonal[i];
    }
    printf("succesfully computed diagonal.\n");
    
    for(int i = 0; i < numrows; i++)
    {
        if(Diagonal[i + firstrow] == 0)
        {
            Diagonal[i + firstrow] = 1;
            baseRows[i][10] = i;
            numElements++;
        }
    }
    int rows[numElements];
    int offsets[numrows];
    offsets[0] = 0;
    for(int i = 0; i < numrows; i++)
    {
        int k = offsets[i];
        for(int l = 0; l < 11; l++)
        {
            if(baseRows[i][l] != -1)
            {
                rows[k] = baseRows[i][l];
                if(rows[k] > N)
                    printf(" error. rows[k] = %i", rows[k]);
                k++;
            }
        }
        if(i != numrows-1)
            offsets[i+1] = k;
    }
    for(int i = 0; i < N; i++)
    {
        if(Diagonal[i] == 0)
            Diagonal[i] = 1;
        Diagonal[i] = 1 / Diagonal[i];
    }

    float u[numrows], res[numrows], tempr[N];
    int tot = 0;

    printf("%i Computed own rows\n", s);
    for(int i = 0; i < numrows; i++)
    {
        int k = rand() % (N*1000);
        u[i] = k;
        tot += k;
    }
    for(int r = 0; r < p; r++)
    {
        if(r != s)
            MPI_Isend(&tot, 1, MPI_INT, r, r, comm, &requests[r]);
    }
    MPI_Barrier(comm);

    for(int r = 0; r < p; r++)
    {
        if(r != s)
        {
            int temptot;
            MPI_Irecv(&temptot, 1, MPI_INT, r, s, comm, &requests[p+r]);
            tot += temptot;
        }
    }
    for(int i = 0; i < numrows; i++)
    {
        u[i] = u[i] / (float)tot;
    }
    printf("%i Computed own u\n", s);
    float prob = 0.85;
    for(int i = 0; i < numrows; i++)
    {
        res[i] = 0;
        tempr[i + firstrow] = u[i] * Diagonal[i + firstrow];
    }
    for(int r = 0; r < p; r++)
        if(r != s) 
        {
            for(int i = 0; i < numrows; i++)
                printf("%i. Outgoing tempr value %f", s, u[i]);
            MPI_Isend(u, numrows, MPI_FLOAT, r, r, comm, &requests[r]);
        }
    MPI_Barrier(comm);
    for(int r = 0; r < p; r++)
    {
        if(r != s) 
        {
            int tempu[numRows(N, p, r)];
            MPI_Irecv(tempu, numRows(N, p, r), MPI_FLOAT, r, s, comm, &requests[p+r]);
            for(int i = 0; i < numRows(N, p, r); i++)
            {
                printf("%i. Incoming tempr value %f\n", s, tempu[i]);
                tempr[i + firstRow(N, p, r)] = tempu[i] * Diagonal[i + firstRow(N, p, r)];
            }
        }
    }
    printf("%i Computed tempr\n", s);
    for(int i = 0; i < N; i++)
        printf("%i . tempr at %i is %f\n", s, i, tempr[i]);
    for(int i = 0; i < numrows; i++)
    {
        int nextOffset;
        if(i == numrows-1)
            nextOffset = numElements;
        else
            nextOffset = offsets[i+1];
        for(int j = offsets[i]; j < nextOffset; j++)
            res[i] += tempr[rows[j]];
        res[i] = res[i] * prob;
        printf("res %i is %f\n", i, res[i]);
        res[i] = 1 - (u[i] - res[i]);
    }
    printf("%i. Computed residual\n", s);
    float norm = 0;
    for(int i = 0; i < numrows; i++)
        norm += res[i]*res[i];

    for(int r = 0; r < p; r++)
        if(r != s)
            MPI_Isend(&norm, 1, MPI_FLOAT, r, r, comm, &requests[r]);
    MPI_Barrier(comm);
    for(int r = 0; r < p; r++)
        if(r != s)
        {
            float temp;
            MPI_Irecv(&temp, 1, MPI_FLOAT, r, s, comm, &requests[p+r]);
            norm += temp;
        }
    norm = sqrt(norm);
    printf("%i. Norm is %f\n", s, norm);

    int t = 0;
    for(int i = 0; i < numrows; i++)
        printf("%i . r %f . u %f\n", firstrow + i, res[i], u[i]);
    while(norm > 0.000001)
    {
        for(int i = 0; i < numrows; i++)
        {
            u[i] += res[i];
            //printf("u = %f\n", u[i]);
        }
        
        //Computed u.
        float newres[numrows];
        for(int i = 0; i  < numrows; i++)
        {
            tempr[i + firstrow] = res[i]*Diagonal[i + firstrow];
            newres[i] = 0;
        }
        for(int r = 0; r < p; r++)
            if(r != s) 
                MPI_Isend(res, numrows, MPI_FLOAT, r, r, comm, &requests[r]);
        MPI_Barrier(comm);
        
        for(int r = 0; r < p; r++)
        {
            if(r != s) 
            {
                int temp[numRows(N, p, r)];
                MPI_Irecv(temp, numRows(N, p, r), MPI_FLOAT, r, s, comm, &requests[p+r]);
                for(int i = 0; i < numRows(N, p, r); i++)
                    tempr[i + firstRow(N, p, r)] = temp[i] * Diagonal[i + firstRow(N, p, r)];;
            }
        }
        MPI_Barrier(comm);
        //Computed tempr.
        for(int i = 0; i < numrows; i++)
        {
            int nextOffset;
            if(i == numrows-1)
                nextOffset = numElements;
            else
                nextOffset = offsets[i+1];
            for(int j = offsets[i]; j < nextOffset; j++)
                newres[i] += tempr[rows[j]];
            res[i] = newres[i] * prob;
        }
        //Computed r.
        for(int i = 0; i < numrows; i++)
            norm += res[i]*res[i];
        for(int r = 0; r < p; r++)
            if(r != s)
                MPI_Isend(&norm, 1, MPI_FLOAT, r, r, comm, &requests[r]);
        MPI_Barrier(comm);
        for(int r = 0; r < p; r++)
            if(r != s)
            {
                float temp;
                MPI_Irecv(&temp, 1, MPI_FLOAT, r, s, comm, &requests[p+r]);
                norm += temp;
            }
        norm = sqrt(norm);
        //printf("%i. Norm in step %i is %f\n", s, t, norm);
        t++;
        if(t > 50)
        {
            break;
        }
    }
}