#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

static bool output = true;

long main(int argc, char **argv)
{
    int p, s;

    /* SPMD part */
    MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD,&p); /* p = number of processors */ 
    MPI_Comm_rank(MPI_COMM_WORLD,&s); /* s = processor number */ 

    int n=1000;
    if (argc>1) n=atoi(argv[1]);

    if(n<0) MPI_Abort(MPI_COMM_WORLD,-1);

    for(long i = 10; i <= n; i *= 10)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if(s == 0)
            printf("Computing n = %i.\n", i);
        computeVector(i, p, s, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
} /* end main */


long numRows(long N, int p, int s)
{
    return (N+p-s-1)/p;
}
long firstRow(long N, int p, int s)
{
    long numrows = numRows(N, p,s);
    long remainder = N % p;
    long firstrow;
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
void computeVector(long N, int p, int s, MPI_Comm comm)
{
    clock_t start, startloop, end;
    start = clock();
    srand(time(0) + s);
    long numrows = (N+p-s-1)/p ;
    long firstrow = firstRow(N, p, s);
    long lastrow = firstrow + numrows - 1;
    long baseRows[numrows][11];
    long numElements = 0;
    MPI_Request requests[2*p];
    for(long i = 0; i < numrows; i++)
    {
        long k = rand() % 10;
        numElements += k + 1;
        for(long l = 0; l < 10; l++)
        {
            if(l <= k)
            {
                baseRows[i][l] = rand() % N;
                if(output)
                    printf(" %ld ", baseRows[i][l]);
            }
            else
            {
                if(output)
                    printf("   ");
                baseRows[i][l] = -1;
            }
        }
        baseRows[i][10] = -1;
        if(output)
            printf("\n");
        if(baseRows[i][10] != -1)
            printf("%i. early invalid value %ld at %ld\n", s, baseRows[i][10], i);
    }
    if(output)
        printf("%i. Generated inlinks\n", s);
    long outgoingDiagonal[p];
    for(int r = 0; r < p; r++)
        outgoingDiagonal[r] = 0;
    for(long i = 0; i < numrows; i++)
        for(int l = 0; l < 10; l++)
            if(baseRows[i][l] != -1)
                outgoingDiagonal[pLoc(N, p, baseRows[i][l])]++;
    if(output)
        printf("%i. Generated local and outgoing diagonals.\n", s);
    
    long outgoingLinks[numElements];
    long outOffsets[p];
    long tempOffsets[p];
    outOffsets[0] = tempOffsets[0] = 0;
    for(int r = 0; r < p - 1; r++)
        outOffsets[r+1] = tempOffsets[r+1] = outOffsets[r] + outgoingDiagonal[r];
    for(long i = 0; i < numrows; i++)
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
    
    long sizes[p];
    for(long r = 0; r < p; r++)
    {
        MPI_Isend(outgoingDiagonal + r, 1, MPI_LONG, r, r, comm, &requests[r]);
        MPI_Irecv(sizes + r, 1, MPI_LONG, r, s, comm, &requests[p+r]);
    }
    MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);
    if(output)
        printf("%i. Sent across sizes. Own size is %ld\n", s, sizes[s]);
    long maxsize = 0;
    for(int r = 0; r < p; r++)
        if(maxsize < sizes[r])
            maxsize = sizes[r];
    long incoming[p*maxsize];
    for(int r = 0; r < p; r++)
        for(int j = 0; j < numElements; j++)
            printf("element %ld\n", outgoingLinks[j]);
    for(long r = 0; r < p; r++)
    {
        MPI_Isend(outgoingLinks + outOffsets[r], outgoingDiagonal[r], MPI_INT, r, r, comm, &requests[r]);
        MPI_Irecv(incoming + r * maxsize, sizes[r] , MPI_LONG, r, s, comm, &requests[p+r]);
    }
    MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);
    if(output)
        printf("%i. Sent across diagonal\n", s);
    long numOutlinks[numrows];
    for(long i = 0; i < numrows; i++)
        numOutlinks[i] = 0;
    for(long r = 0; r < p; r++)
        for(long i = 0; i < sizes[r]; i++)
        {
            if(incoming[i + maxsize * r] > N)
                printf("%i. Incoming value %ld", s, incoming[i + maxsize * r]);
            numOutlinks[incoming[i + maxsize * r] - firstrow]++;
        }
    /*
    for(long i = 0; i < N; i++)
        localDiagonal[i] = 0;
    for(long i = 0; i < numrows; i++)
    {
        for(long l = 0; l < 11; l++)
            if(baseRows[i][l] != -1)
                localDiagonal[baseRows[i][l]]++;
        if(baseRows[i][10] != -1)
            printf("%i. middle invalid value %ld at %ld\n", s, baseRows[i][10], i);
    }
    int temp[p*N];
    
    
    for(long r = 0; r < p; r++)
    {
        MPI_Isend(localDiagonal, N, MPI_INT, r, r, comm, &requests[r]);
        MPI_Irecv(temp + r*N, N, MPI_INT, r, s, comm, &requests[p+r]);
    }
    MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);
    
    long numOutlinks[N];
    for(int i = 0; i < N; i++)
    {
        numOutlinks[i] = 0;
        
        for(int r = 0; r < p; r++)
            numOutlinks[i] += temp[N*r+i];
    }
    */
    if(output)
        printf("%i. Generated outlinks.\n", s);
    
    long selfLinks = 0;
    for(long i = 0; i < numrows; i++)
    {
        if(baseRows[i][10] != -1)
            printf("%i. invalid value %ld at %ld\n", s, baseRows[i][10], i);
        if(numOutlinks[i] == 0)
        {
            selfLinks++;
            numElements++;
            baseRows[i][10] = i + firstrow;
            numOutlinks[i] = 1;
        }
    }
    if(output)
        printf("%i. Added additional selflinks.\n", s);
    long rows[numElements];
    long offsets[numrows];
    offsets[0] = 0;
    for(long i = 0; i < numrows; i++)
    {
        long k = offsets[i];
        for(int l = 0; l < 11; l++)
        {
            if(baseRows[i][l] != -1)
            {
                rows[k] = baseRows[i][l];
                if(rows[k] > N)
                    printf("%i. Invalid row value %ld at %ld\n", s, rows[k], k);
                k++;
            }
        }
        if(i + 1 != numrows)
            offsets[i+1] = k;
        else
            if(k != numElements)
                printf("%i. k should be %ld, but is %ld\n Selflinks = %ld\n", s, numElements, k, selfLinks);
    }
    float Diagonal[numrows];
    for(long i = 0; i < numrows; i++)
    {
        if(numOutlinks[i] == 0)
            numOutlinks[i] = 1;
        Diagonal[i] = 1 / (float)numOutlinks[i];
        if(output)
            printf("%i. Diagonal at %ld is %f\n", s, i, Diagonal[i]);
    }
    if(output)
        printf("Computed stochastic row Matrix.\n");
    float u[numrows], res[numrows], tempr[N];
    int tot = rand() % 1000;
    long loctot = 1;
    
    for(long i = 0; i < numrows; i++)
    {
        int k = rand() % 1000;
        u[i] = k * (float)tot;
        loctot += k;
    }
    for(long i = 0; i < numElements; i++)
        if(rows[i] > N)
            printf("%i. OLD Rows[%ld] value is %ld\n", s, i, rows[i]);
    
    MPI_Barrier(comm);
    
    int tmp[p];
    for(long r = 0; r < p; r++)
    {
        MPI_Isend(&tot, 1, MPI_INT, r, r, comm, &requests[r]);
        MPI_Irecv(&tmp, 1, MPI_INT, r, s, comm, &requests[p+r]);
    }
    MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);
    tot = 0;
    for(int r = 0; r < p; r++)
        tot += tmp[r];
    for(long i = 0; i < numElements; i++)
        if(rows[i] > N)
            printf("%i. Rows[%ld] value is %ld\n", s, i, rows[i]);
    
    for(long i = 0; i < numrows; i++)
        u[i] = u[i] / (float)tot /(float)loctot;
    if(output)
        printf("%i Computed own u\n", s);
    long t = 0;
    float norms[1000];
    for(long i = 0; i < 1000; i++)
        norms[i] = -1;
    float prob = 0.85;
    float newres[numrows];

    for(long i = 0; i < numrows; i++)
    {
        res[i] = 0;
        newres[i] = u[i] * Diagonal[i];
    }
    float tempu[p * numRows(N,p,0)];
    for(long r = 0; r < p; r++)
    {
        MPI_Isend(newres, numrows, MPI_FLOAT, r, r, comm, &requests[r]);
        MPI_Irecv(tempu + r * numRows(N,p,0), numRows(N, p, r), MPI_FLOAT, r, s, comm, &requests[p+r]);
    }
    MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);
    for(int r = 0; r < p; r++)
        for(long i = 0; i < numRows(N, p, r); i++)
            tempr[i + firstRow(N, p, r)] = tempu[i + r * numRows(N, p, 0)]; 

    if(output)
        printf("%i Computed tempr\n", s);
    for(long i = 0; i < numrows; i++)
    {
        long nextOffset;
        if(i == numrows-1)
            nextOffset = numElements;
        else
            nextOffset = offsets[i+1];
        for(long j = offsets[i]; j < nextOffset; j++) 
        {
            if(j > numElements)
                printf("%i. Invalid j = %ld\n", s, j);
            if(rows[j] > N)
            {
                printf("%i. Invalid rows value %ld at position %ld\n", s, rows[j], j);
                printf("%i. N = %ld\n", s, N);
            }
            res[i] += tempr[rows[j]];
        }
        res[i] = res[i] * prob;
        res[i] = 1 - (u[i] - res[i]);
    }
    if(output)
        printf("%i. Computed initial residual\n", s);
    
    float norm = 0;
    for(long i = 0; i < numrows; i++)
        norm += res[i]*res[i];
    float locnorms[p];
    for(long r = 0; r < p; r++)
    {
        MPI_Isend(&norm, 1, MPI_FLOAT, r, r, comm, &requests[r]);
        MPI_Irecv(locnorms + r, 1, MPI_FLOAT, r, s, comm, &requests[p+r]);
    }
    MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);
    norm = 0;
    for(int r = 0; r < p; r++)
        norm += locnorms[r];
    
    norm = sqrt(norm);
    if(output)
        printf("%i. Norm is %f\n", s, norm);
    startloop = clock();
    
    while(norm > 0.000001)
    {
        for(long i = 0; i < numrows; i++)
            u[i] += res[i];

        //Computed u.
        for(long i = 0; i  < numrows; i++)
        {
            res[i] = res[i] * Diagonal[i];
            newres[i] = 0;
        }

        for(long r = 0; r < p; r++)
        {
            MPI_Isend(res, numrows, MPI_FLOAT, r, r, comm, &requests[r]);
            MPI_Irecv(tempu + r * numRows(N, p, 0) , numRows(N, p, r), MPI_FLOAT, r, s, comm, &requests[p + r]);
        }
        MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);
        
        for(int r = 0; r < p; r++)
        {
            for(long i = 0; i < numRows(N,p,r); i++)
                    tempr[i + firstRow(N, p, r)] = tempu[i + r * numRows(N,p,r)];
        }
        //Computed tempr.
        for(long i = 0; i < numrows; i++)
        {
            long nextOffset;
            if(i == numrows-1)
                nextOffset = numElements;
            else
                nextOffset = offsets[i+1];
            for(long j = offsets[i]; j < nextOffset; j++)
                newres[i] += (float)tempr[rows[j]];
            res[i] = newres[i] * prob;
        }
        //Computed r.
        norm = 0;
        for(long i = 0; i < numrows; i++)
            norm += res[i]*res[i];
        for(long r = 0; r < p; r++)
        {
            MPI_Isend(&norm, 1, MPI_FLOAT, r, r, comm, &requests[r]);
            MPI_Irecv(locnorms + r, 1, MPI_FLOAT, r, s, comm, &requests[p + r]);
        }
        MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);
        norm = 0;
        for(long r = 0; r < p; r++)
            norm += locnorms[r];
        norm = sqrt(norm);
        norms[t] = norm;
        if(output && t < 5)
            printf("%i. Norm in step %i is %f\n", s, t, norm);
        t++;
        if(t > 1000)
        {
            printf("%i. Loop break at t = %i. Norm is %f\n", s, t, norm);
            break;
        }
    }
    
    end = clock();
    if(s == 0)
    {
        float tottime = ((float)(end - start)) / CLOCKS_PER_SEC;
        float initialtime = ((float)(startloop - start)) / CLOCKS_PER_SEC;
        float looptime = ((float)(end - startloop)) / CLOCKS_PER_SEC;
        printf("total time is %f\n", tottime);
        printf("initial time is %f\n", initialtime);
        printf("loop time is %f\n", looptime);
    }
}