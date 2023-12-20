#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

static bool output = true;

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

    for(int i = 10; i <= n; i *= 10)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if(s == 0)
            printf("Computing n = %i.\n", i);
        computeVector(i, p, s, MPI_COMM_WORLD);
    }

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

int pLoc(int N, int p, int j)
{
    int remainder = N % p;
    int numrows = N / p;
    if(j < remainder * (numrows + 1))
        return j / (numrows + 1);
    else
        return j / numrows - remainder;
}

//n is number of rows/columns.
//p is number of processors.
//s is own rank.
void computeVector(int N, int p, int s, MPI_Comm comm)
{
    clock_t start, startloop, end;
    start = clock();
    srand(time(0) + s);
    int numrows = (N+p-s-1)/p ;
    int firstrow = firstRow(N, p, s);
    int lastrow = firstrow + numrows - 1;
    int baseRows[numrows][10];
    for(int i = 0; i < numrows; i++)
    {
        int k = rand() % 10;
        for(int l = 0; l < 10; l++)
        {
            if(l <= k)
            {
                baseRows[i][l] = rand() % N;
                if(output && i < 2)
                    printf(" %i ", baseRows[i][l]);
            }
            else
            {
                if(output && i < 2)
                    printf("   ");
                baseRows[i][l] = -1;
            }
        }
        if(output && i < 2)
            printf("\n");
    }
    if(output)
        printf("%i. Generated inlinks\n", s);
    int localDiagonal[numrows];
    int outgoingDiagonal[p];
    int numElements = 0;
    for(int i = 0; i < N; i++)
        localDiagonal[i] = 0;
    for(int r = 0; r < p; r++)
        outgoingDiagonal[r] = 0;
    for(int i = 0; i < numrows; i++)
        for(int l = 0; l < 10; l++)
            if(baseRows[i][l] != -1)
            {
                numElements++;
                int j = baseRows[i][l];
                if(j >= firstrow && j <= lastrow)
                    localDiagonal[j - firstrow]++;
                else
                    outgoingDiagonal[pLoc(N, p, j)]++;
            }
    if(output)
        printf("%i. Generated local and outgoing diagonals.\n", s);
    int outgoingLinks[numElements];
    int outOffsets[p];
    int tempOffsets[p];
    outOffsets[0] = tempOffsets[0] = 0;
    for(int r = 0; r < p - 1; r++)
        outOffsets[r+1] = tempOffsets[r+1] = outOffsets[r] + outgoingDiagonal[r];
    for(int i = 0; i < numrows; i++)
        for(int l = 0; l < 10; l++)
            if(baseRows[i][l] != -1)
            {
                int j = baseRows[i][l];
                if(j < firstrow || j > lastrow)
                {
                    outgoingLinks[tempOffsets[pLoc(N, p, j)]] = j;
                    tempOffsets[pLoc(N, p, j)]++;
                }
            }
    if(output)
        printf("%i Generated Outgoing links.\n", s);
    MPI_Request requests[2*p];
    for(int r = 0; r < p; r++)
    {
        if(r == s) continue;
        MPI_Isend(&outgoingDiagonal[r], 1, MPI_INT, r, r, comm, &requests[r]);
        MPI_Isend(&outgoingLinks[outOffsets[r]], outgoingDiagonal[r], MPI_INT, r, r, comm, &requests[r]);
    }
    MPI_Barrier(comm);
    for(int r = 0; r < p; r++)
    {
        if(r == s) continue;
        int size;
        
        MPI_Irecv(&size, 1, MPI_INT, r, s, comm, &requests[p+r]);
        int incoming[size];
        MPI_Irecv(incoming, size, MPI_INT, r, s, comm, &requests[p+r]);
        for(int i = 0; i < size; i++)
            localDiagonal[incoming[i] - firstrow]++;
    }
    MPI_Barrier(comm);
    if(output)
        printf("%i. Generated outlinks.\n", s);
    
    if(output)
        printf("%i. Added additional selflinks.\n", s);
    int rows[numElements];
    int offsets[numrows];
    offsets[0] = 0;
    for(int i = 0; i < numrows; i++)
    {
        int k = offsets[i];
        for(int l = 0; l < 10; l++)
        {
            if(baseRows[i][l] != -1)
            {
                rows[k] = baseRows[i][l];
                k++;
            }
        }
        if(i != numrows-1)
            offsets[i+1] = k;
    }
    
    float Diagonal[numrows];
    for(int i = 0; i < numrows; i++)
    {
        if(localDiagonal[i] == 0)
            localDiagonal[i] = 1;
        Diagonal[i] = 1 / (float)localDiagonal[i];
        if(output && i < 10)
            printf("%i. Diagonal at %i is %f\n", s, i + firstrow, Diagonal[i]);
    }
    if(output)
        printf("Computed stochastic row Matrix.\n");
    float u[numrows], res[numrows], tempr[N];
    int tot = 0;

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
        u[i] /=  (float)tot;
    if(output)
        printf("%i Computed own u\n", s);
    int t = 0;
    float norms[1000];
    for(int i = 0; i < 1000; i++)
        norms[i] = -1;
    float prob = 0.85;

    
    for(int i = 0; i < numrows; i++)
    {
        res[i] = 0;
        u[i] = u[i] * Diagonal[i];
        tempr[firstrow + i] = u[i];
    }
    for(int r = 0; r < p; r++)
        if(r != s) 
            MPI_Isend(u, numrows, MPI_FLOAT, r, r, comm, &requests[r]);
    MPI_Barrier(comm);
    for(int r = 0; r < p; r++)
        if(r != s) 
            MPI_Irecv(&tempr[firstRow(N, p, r)], numRows(N, p, r), MPI_FLOAT, r, s, comm, &requests[p+r]);
    if(output)
        printf("%i Computed tempr\n", s);
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
        res[i] = 1 - (u[i] - res[i]);
    }
    if(output)
        printf("%i. Computed initial residual\n", s);
    
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
    if(output)
        printf("%i. Norm is %f\n", s, norm);
    startloop = clock();
    while(norm > 0.000001)
    {
        for(int i = 0; i < numrows; i++)
            u[i] += res[i];

        //Computed u.
        float newres[numrows];
        for(int i = 0; i  < numrows; i++)
        {
            res[i] *= Diagonal[i];
            tempr[i + firstrow] = res[i];
            newres[i] = 0;
        }
        for(int r = 0; r < p; r++)
            if(r != s) 
                MPI_Isend(res, numrows, MPI_FLOAT, r, r, comm, &requests[r]);
        MPI_Barrier(comm);
        
        for(int r = 0; r < p; r++)
            if(r != s) 
                MPI_Irecv(&tempr[firstRow(N, p, r)], numRows(N, p, r), MPI_FLOAT, r, s, comm, &requests[p+r]);

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
                newres[i] += (float)tempr[rows[j]];
            res[i] = newres[i] * prob;
        }
        //Computed r.
        norm = 0;
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
        norms[t] = norm;
        if(output && t < 1)
            printf("%i. Norm in step %i is %f\n", s, t, norm);
        t++;
        if(t > 1000)
        {
            printf("broke at norm size %f", norm);
            break;
        }
    }
    end = clock();
    if(s == 0)
    {
        printf("Finished in %i steps.\n", t);
        float tottime = ((float)(end - start)) / CLOCKS_PER_SEC;
        float initialtime = ((float)(startloop - start)) / CLOCKS_PER_SEC;
        float looptime = ((float)(end - startloop)) / CLOCKS_PER_SEC;
        printf("total time is %f\n", tottime);
        printf("initial time is %f\n", initialtime);
        printf("loop time is %f\n", looptime);
    }
}