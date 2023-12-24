#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>


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

    //Start iterating at $N = 10. Repeatedly keep increasing N by a factor 10 until the maximum is reached.
    for(long i = 10; i <= n; i *= 10)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if(s == 0)
            printf("Computing n = %i.\n", i);
        if(i * 10 > n)
            i = n;
        computeVector(i, p, s, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
} /* end main */

//Gives the number of rows on processor s with a total of N rows and p processors. 
long numRows(long N, int p, int s)
{
    return (N + p - s - 1)/p;
}
//Gives the global index of the first row of processor s when there are a total of N rows and p processors.
long firstRow(long N, int p, int s)
{
    long numrows = numRows(N, p,s);
    long remainder = N % p;
    long firstrow;
    if( s < remainder) 
        //numrows is the same as all earlier processors.
        firstrow = s * numrows;
    else
        //The first remainder processors have 1 more element.
        firstrow = s * numrows + remainder;
    return firstrow;
}

//Gives the processor number that stores the row with global index j when there are N rows and p processors. 
int pLoc(int N, int p, int j)
{
    int remainder = N % p;
    int numrows = N / p;
    if(j < remainder * (numrows + 1))
        return j / (numrows + 1);
    else
        return  (j - remainder) / numrows;
        //Equivalent to: remainder + (j - remainder * (numrows + 1)/numrows
}

//n is number of rows/columns.
//p is number of processors.
//s is own rank.
void computeVector(long N, int p, int s, MPI_Comm comm)
{
    //Clocks for measuring the time
    clock_t start, startloop, end;
    start = clock();
    //Set random seed
    srand(time(0) + s);
    long numrows = (N+p-s-1)/p ;
    long firstrow = firstRow(N, p, s);
    //Block partition, so last row is determined from firstrow and numrows.
    long lastrow = firstrow + numrows - 1;
    long baseRows[numrows][11];
    long numElements = 0;
    float prob = 0.85;
    MPI_Request requests[2*p];
    //Creates baseRows. Here positions with -1 will be removed later in the formation of rows.
    for(long i = 0; i < numrows; i++)
    {
        long k = rand() % 10;
        numElements += k + 1;
        for(long l = 0; l < 11; l++)
        {
            if(l <= k)
                baseRows[i][l] = rand() % N;
            else 
                baseRows[i][l] = -1;
        }
    }
    //Keep track of how many elements will be sent to each other processor.
    long outgoingDiagonal[p];
    for(int r = 0; r < p; r++)
        outgoingDiagonal[r] = 0;
    for(long i = 0; i < numrows; i++)
        for(int l = 0; l < 10; l++)
            if(baseRows[i][l] != -1)
                outgoingDiagonal[pLoc(N, p, baseRows[i][l])]++;

    //Keep track which elements will be sent to each other processor. Store this in CRS format.
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
                int r = pLoc(N, p, baseRows[i][l]);
                outgoingLinks[tempOffsets[r]] = baseRows[i][l];
                tempOffsets[r]++;
            }
    
    //Send a message indicating the number of elements that each processor will receive. Store the result in sizes.
    long sizes[p];
    for(long r = 0; r < p; r++)
    {
        MPI_Isend(outgoingDiagonal + r, 1, MPI_LONG, r, r, comm, &requests[r]);
        MPI_Irecv(sizes + r, 1, MPI_LONG, r, s, comm, &requests[p+r]);
    }
    MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);long maxsize = 0;
    //Find the largest number of elements we will be receiving.
    for(int r = 0; r < p; r++)
        if(maxsize < sizes[r])
            maxsize = sizes[r];
    long incoming[p*maxsize];
    //Get elements from the other processors. Some spots in incoming may stay empty.
    for(long r = 0; r < p; r++)
    {
        MPI_Isend(outgoingLinks + outOffsets[r], outgoingDiagonal[r], MPI_LONG, r, r, comm, &requests[r]);
        MPI_Irecv(incoming + r * maxsize, sizes[r] , MPI_LONG, r, s, comm, &requests[p+r]);
    }
    MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);

    //Compute diagonal.
    long numOutlinks[numrows];
    for(long i = 0; i < numrows; i++)
        numOutlinks[i] = 0;
    for(long r = 0; r < p; r++)
        for(long i = 0; i < sizes[r]; i++)
            numOutlinks[incoming[i + maxsize * r] - firstrow]++;

    
    //Add self-links in baseRows and the diagonal.
    long selfLinks = 0;
    for(long i = 0; i < numrows; i++)
        if(numOutlinks[i] == 0)
        {
            selfLinks++;
            numElements++;
            baseRows[i][10] = i + firstrow;
            numOutlinks[i] = 1;
        }
    //Prune baseRows and store result in CRS format.
    long rows[numElements];
    long offsets[numrows];
    offsets[0] = 0;
    for(long i = 0; i < numrows; i++)
    {
        long k = offsets[i];
        for(int l = 0; l < 11; l++)
            if(baseRows[i][l] != -1)
            {
                rows[k] = baseRows[i][l];
                k++;
            }
        if(i + 1 != numrows)
            offsets[i+1] = k;
    }
    //Determine inverse diagonal.
    float Diagonal[numrows];
    for(long i = 0; i < numrows; i++)
        Diagonal[i] = 1 / (float)numOutlinks[i];


    float u[numrows], res[numrows], tempr[N];
    long loctot = 0;
    long globtot;
    
    //Give each u[i] a random value.
    for(long i = 0; i < numrows; i++)
    {
        int k = rand() % 1000;
        u[i] = k;
        loctot += k;
    }
    MPI_Allreduce(&loctot, &globtot, 1, MPI_LONG, MPI_SUM, comm);
    //Normalize u[i]
    for(long i = 0; i < numrows; i++)
        u[i] = u[i] / (float)globtot;

    
    //Determine initial residual.
    float newres[numrows];
    for(long i = 0; i < numrows; i++)
    {
        res[i] = 0;
        newres[i] = u[i] * Diagonal[i];
    }
    float tempu[p * numRows(N,p,0)];
    //Send newres to each other processor.
    for(long r = 0; r < p; r++)
    {
        MPI_Isend(newres, numrows, MPI_FLOAT, r, r, comm, &requests[r]);
        MPI_Irecv(tempu + r * numRows(N, p, 0), numRows(N, p, r), MPI_FLOAT, r, s, comm, &requests[p+r]);
    }
    MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);
    for(int r = 0; r < p; r++)
        for(long i = 0; i < numRows(N, p, r); i++)
            tempr[i + firstRow(N, p, r)] = tempu[i + r * numRows(N, p, 0)]; 

    //Compute initial residual.
    for(long i = 0; i < numrows; i++)
    {
        long nextOffset;
        if(i == numrows-1)
            nextOffset = numElements;
        else
            nextOffset = offsets[i+1];
        for(long j = offsets[i]; j < nextOffset; j++) 
            res[i] += tempr[rows[j]];

        res[i] = res[i] * prob;
        res[i] = 1 + res[i] - u[i];
    }
    
    //Compute initial norm.
    float locnorm = 0;
    float globnorm;
    for(long i = 0; i < numrows; i++)
        locnorm += res[i]*res[i];

    MPI_Allreduce(&locnorm, &globnorm, 1, MPI_FLOAT, MPI_SUM, comm);

    
    float norm = sqrt(globnorm);
    long t = 0;

    //Start loop
    startloop = clock();
    while(norm > 0.000001)
    {
        //compute new u
        for(long i = 0; i < numrows; i++)
            u[i] += res[i];

        //Computed new r * D^-1.
        for(long i = 0; i  < numrows; i++)
        {
            res[i] = res[i] * Diagonal[i];
            newres[i] = 0;
        }
        //Communicate r * D^-1 to all processors 
        for(long r = 0; r < p; r++)
        {
            MPI_Isend(res, numrows, MPI_FLOAT, r, r, comm, &requests[r]);
            MPI_Irecv(tempu + r * numRows(N, p, 0) , numRows(N, p, r), MPI_FLOAT, r, s, comm, &requests[p + r]);
        }
        MPI_Waitall(2*p, requests,MPI_STATUSES_IGNORE);
        
        for(int r = 0; r < p; r++)
        {
            for(long i = 0; i < numRows(N,p,r); i++)
                tempr[i + firstRow(N, p, r)] = tempu[i + r * numRows(N, p, 0)];
        }
        //Compute next r.
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
        //Compute new norm.
        norm = 0;
        locnorm = 0;
        for(long i = 0; i < numrows; i++)
            locnorm += res[i]*res[i];
        MPI_Allreduce(&locnorm, &globnorm, 1, MPI_FLOAT, MPI_SUM, comm);

        norm = sqrt(globnorm);
        t++;
        //backup loop break to ensure the program terminates.
        if(t > 1000)
        {
            printf("%i. Loop break at t = %i. Norm is %f\n", s, t, norm);
            break;
        }
    }
    //Final time measurements
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