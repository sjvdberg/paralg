#include <stdio.h>
#include <stdbool.h>
#include <math.h>

static long N = 10; 

void main()
{
    printf("N = %i\n", N);
    int baseRows[N][11];
    for (int i = 0; i < N; i++)
    {
        int k = rand() % 10;
        for(int l = 0; l < 11; l++)
        {
            if(l <= k)
            {
                baseRows[i][l] = rand() % N;
            }
            else
            {
                baseRows[i][l] = -1;
            }
        }
    }

    printf("Generated inlinks.\n");

    //Get number of Outlinks.
    int numOutlinks[N];
    for(int i = 0; i < N; i++)
        numOutlinks[i] = 0;
    for(int i = 0; i < N; i++)
        for(int j = 0; j < 11; j++)
            if(baseRows[i][j] != -1)
                numOutlinks[baseRows[i][j]]++;

    printf("Computed outlinks.\n");

    int numElements = 0;
    for (int i = 0; i < N; i++)
    {
        if(numOutlinks[i] == 0)
        {
            numElements++;
            baseRows[i][10] = i;
            numOutlinks[i] = 1;
        }
        else
        {
            numElements += numOutlinks[i];
        }
    }

    printf("Added additional outlinks.\n");
    int rows[numElements];
    int offsets[N];
    offsets[0] = 0;

    for (int i = 0; i < N; i++)
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
        if(i != N - 1) 
            offsets[i+1] = k;
    }

    float diagonal[N];
    for(int i = 0; i < N; i++)
        diagonal[i] = 1 / (float)numOutlinks[i];

    printf("Computed stochastic row Matrix.\n");
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

    printf("Computed initial u.\n");
    int t = 0;
    float norms[1000];
    for(int i = 0; i <1000; i++)
        norms[i] = -1;
    float p = 0.85;

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
        r[i] = 1 - (u[i] - r[i]);
    }
    printf("Computed initial r.\n");

    float norm = 0;
    for(int i = 0; i < N; i++)
        norm += r[i]*r[i];
    norm = sqrt(norm);
    for(int i = 0; i < N; i++)
        printf("u %f . r %f\n", u[i], r[i]);
    while(norm > 0.000001)
    {
        for(int i = 0; i < N; i++)
            u[i] += r[i];
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
        norm = 0;
        for(int i = 0; i< N; i++)
            norm += r[i]*r[i];
        norm = sqrt(norm);
        norms[t] = norm;
        printf("Computed u and r in step %i, current norm is %f\n", t, norm);
        t++;
    }
    float totNormChanges = 0;
    for(int i = 0; i < t -1; i++)
    {
        totNormChanges += (norms[i+1]/norms[i]);
    }
    totNormChanges = totNormChanges / t;
    printf("average norm change is %f\n", totNormChanges);
}
