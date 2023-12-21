#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

static bool output = false;

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
    int baseRows[numrows][11];
    int numElements = 0;
    for(int i = 0; i < numrows; i++)
    {
        int k = rand() % 10;
        numElements += k + 1;
        for(int l = 0; l < 11; l++)
        {
            if(l <= k)
            {
                baseRows[i][l] = rand() % N;
                if(output)
                    printf(" %i ", baseRows[i][l]);
            }
            else
            {
                if(output)
                    printf("   ");
                baseRows[i][l] = -1;
            }
        }
        if(output)
            printf("\n");
    }
    if(output)
        printf("%i. Generated inlinks\n", s);
    int localDiagonal[N];
    for(int i = 0; i < N; i++)
        localDiagonal[i] = 0;
    for(int i = 0; i < numrows; i++)
        for(int l = 0; l < 11; l++)
            if(baseRows[i][l] != -1)
            {
                localDiagonal[baseRows[i][l]]++;
            }
    int numOutlinks[N];
    MPI_Request requests[2*p];
    for(int r = 0; r < p; r++)
    {
        if(r == s)
            for(int i = 0; i < N; i++)
                numOutlinks[i] = localDiagonal[i];
        else
            MPI_Isend(localDiagonal, N, MPI_INT, r, r, comm, &requests[r]);
    }
    MPI_Barrier(comm);
    for(int r = 0; r < p; r++)
    {
        if(r == s) continue;
        MPI_Irecv(localDiagonal, N, MPI_INT, r, s, comm, &requests[p+r]);
        for(int i = 0; i < N; i++)
            numOutlinks[i] += localDiagonal[i];
    }
    if(output)
        printf("%i. Generated outlinks.\n", s);
    
    int selfLinks = 0;
    for(int i = 0; i < numrows; i++)
    {
        if(numOutlinks[i + firstrow] == 0)
        {
            selfLinks++;
            numElements++;
            baseRows[i][10] = i + firstrow;
            numOutlinks[i + firstrow] = 1;
        }
    }
    if(output)
        printf("%i. Added additional selflinks.\n", s);
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
                k++;
            }
        }
        if(i != numrows-1)
            offsets[i+1] = k;
        else
            if(k != numElements)
                printf("%i. k should be %i, but is %i\n Selflinks = %i\n", s, numElements, k, selfLinks);
    }
    
    float Diagonal[N];
    for(int i = 0; i < N; i++)
    {
        if(numOutlinks[i] == 0)
            numOutlinks[i] = 1;
        Diagonal[i] = 1 / (float)numOutlinks[i];
        if(output)
            printf("%i. Diagonal at %i is %f\n", s, i, Diagonal[i]);
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
        tempr[i + firstrow] = u[i] * Diagonal[i + firstrow];
    }
    for(int r = 0; r < p; r++)
        if(r != s) 
            MPI_Isend(u, numrows, MPI_FLOAT, r, r, comm, &requests[r]);
    MPI_Barrier(comm);
    for(int r = 0; r < p; r++)
    {
        if(r != s) 
        {
            float tempu[numRows(N, p, r)];
            MPI_Irecv(tempu, numRows(N, p, r), MPI_FLOAT, r, s, comm, &requests[p+r]);
            for(int i = 0; i < numRows(N, p, r); i++)
                tempr[i + firstRow(N, p, r)] = tempu[i] * Diagonal[i + firstRow(N, p, r)];
        }
    }
    
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
        {
            if(j > numElements)
                printf("%i. Invalid j = %i\n", s, j);
            if(rows[j] > N)
                printf("%i. Invalid j value %i at position %i\n", s, rows[j], j);
            res[i] += tempr[rows[j]];
        }
        res[i] = res[i] * prob;
        res[i] = 1 - (u[i] - res[i]);
    }
    if(output)
        printf("%i. Computed initial residual\n", s);
    /*
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
            tempr[i + firstrow] = res[i] * Diagonal[i+firstrow];
            newres[i] = 0;
        }
        for(int r = 0; r < p; r++)
            if(r != s) 
                MPI_Isend(res, numrows, MPI_FLOAT, r, r, comm, &requests[r]);
        MPI_Barrier(comm);
        
        for(int r = 0; r < p; r++)
            if(r != s) 
            {
                float temp[numRows(N,p,r)];
                MPI_Irecv(temp, numRows(N, p, r), MPI_FLOAT, r, s, comm, &requests[p+r]);
                for(int i = 0; i < numRows(N,p,r); i++)
                    tempr[i + firstRow(N, p, r)] = temp[i] * Diagonal[i + firstRow(N, p, r)];
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
        if(output)
            printf("%i. Norm in step %i is %f\n", s, t, norm);
        t++;
        if(t > 1000)
        {
            printf("%i. Loop break at t = %i. Norm is %f\n", s, t, norm);
        }
    }
    */
    end = clock();
    if(s == 0)
    {
        float tottime = ((float)(end - start)) / CLOCKS_PER_SEC;
        //float initialtime = ((float)(startloop - start)) / CLOCKS_PER_SEC;
        //float looptime = ((float)(end - startloop)) / CLOCKS_PER_SEC;
        printf("total time is %f\n", tottime);
        //printf("initial time is %f\n", initialtime);
        //printf("loop time is %f\n", looptime);
    }
}