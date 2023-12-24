#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

static bool output = false;

void main()
{
    //iterate through increasing values of N until a certain limit is reached.
    int N = 10000000;
    for(int i = 10; i <= N; i *= 10)
        Serial(i);
}

void Serial(long N)
{
    //Use clocks for time measurement.
    clock_t start, startloop, end;
    start = clock();
    float p = 0.85;
    //Counter for the number of elements in the graph.
    int numElements = 0;
    printf("N = %i\n", N);
    //Determine initial rows. -1 values will be removed later.
    int baseRows[N][11];
    for (int i = 0; i < N; i++)
    {
        int k = rand() % 10;
        numElements += k+1;
        for(int l = 0; l < 11; l++)
        {
            if(l <= k)
                baseRows[i][l] = rand() % N;
            else
                baseRows[i][l] = -1;
        }
    }

    //Get number of Outlinks from each page.
    int numOutlinks[N];
    for(int i = 0; i < N; i++)
        numOutlinks[i] = 0;
    for(int i = 0; i < N; i++)
        for(int j = 0; j < 11; j++)
            if(baseRows[i][j] != -1)
                numOutlinks[baseRows[i][j]]++;

    //Add self-links.
    for (int i = 0; i < N; i++)
        if(numOutlinks[i] == 0)
        {
            numElements++;
            baseRows[i][10] = i;
            numOutlinks[i] = 1;
        }

    //Prune baseRows and store result in CRS format.
    int rows[numElements];
    int offsets[N];
    offsets[0] = 0;
    for (int i = 0; i < N; i++)
    {
        int k = offsets[i];
        for(int l = 0; l < 11; l++)
            if(baseRows[i][l] != -1)
            {
                rows[k] = baseRows[i][l];
                k++;
            }
        if(i != N - 1) 
            offsets[i+1] = k;
    }
    //Compute inverse of diagonal.
    float diagonal[N];
    for(int i = 0; i < N; i++)
        diagonal[i] = 1 / (float)numOutlinks[i];


    //Compute initial u.
    float u[N], r[N], tempr[N];
    int tot = 0;
    for(int i = 0; i < N; i++)
    {
        int k = rand() % N*1000;
        u[i] = k;
        tot += k;
    }
    for(int i = 0; i < N; i++)
        u[i] /= tot;

    //Compute initial r.
    for(int i = 0; i < N; i++)
    {
        r[i] = 0;
        tempr[i] = u[i] * diagonal[i];
    }
    for(int i = 0; i < N; i++)
    {
        int nextOffset;
        if(i == N-1)
            nextOffset = numElements;
        else
            nextOffset = offsets[i+1];
        for(int j = offsets[i]; j < nextOffset; j++)
            r[i] += tempr[rows[j]];
        r[i] = r[i] * p;
        r[i] = 1 + r[i] - u[i];
    }

    //Compute intial norm
    float norm = 0;
    for(int i = 0; i < N; i++)
        norm += r[i]*r[i];
    norm = sqrt(norm);

    //Create table to keep track of norms in each iteration.
    int norms[1000];
    for(int i = 0; i < 1000; i++)
        norms[i] = -1;
    norms[0] = norm;
    startloop = clock();
    int t = 0;
    while(norm > 0.000001)
    {
        //Compute new u.
        for(int i = 0; i < N; i++)
            u[i] += r[i];
        //Compute new r.
        for(int i = 0; i < N; i++)
        {
            tempr[i] = r[i]*diagonal[i];
            r[i] = 0;
        }
        for(int i = 0; i < N; i++)
        {
            int nextOffset;
            if(i == N-1)
                nextOffset = numElements;
            else
                nextOffset = offsets[i+1];
            for(int j = offsets[i]; j < nextOffset; j++)
                r[i] +=(float)tempr[rows[j]];
            r[i] = r[i] * p;
        }
        //Compute new norm
        norm = 0;
        for(int i = 0; i< N; i++)
            norm += r[i]*r[i];
        norm = sqrt(norm);
        t++;
        norms[t] = norm;
    }
    //compute average norm changes
    float totNormChanges = 0;
    for(int i = 0; i < t -1; i++)
        totNormChanges += (norms[i+1]/norms[i]);
    totNormChanges = totNormChanges / t;
    printf("average norm change is %f\n", totNormChanges);

    //Do time measurements.
    end = clock();
    float tottime = ((float)(end - start)) / CLOCKS_PER_SEC;
    float initialtime = ((float)(startloop - start)) / CLOCKS_PER_SEC;
    float looptime = ((float)(end - startloop)) / CLOCKS_PER_SEC;
    printf("total time is %f\n", tottime);
    printf("initial time is %f\n", initialtime);
    printf("loop time is %f\n", looptime);
}
