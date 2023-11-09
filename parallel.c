#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <bsp.h>

static long P; // number of processors requested
static long N; // upper bound on the primes.
static long Num; // number of elements per processor

void sieve()
{
    bsp_begin(P);
    long p = P;
    long s = bsp_pid();
    long n = N;

    int numvalues, startvalue;
    numvalues = Num;
    startvalue = s * numvalues;

    bool values[Num];
    for(int i = 0; i < numvalues; i++)
        values[i] = false;
    int lowestindex = 0;
    int currentprime = 2;
    if(s == 0)
    {
        values[0] = true;
        values[1] = true;
        lowestindex = 2;
    }
    double *Lowest= malloc(n*sizeof(long));
    bsp_push_reg(Lowest,p*sizeof(long));
    bsp_sync();
    printf("start sieve with %d processors on %d.\n", p, s);


    while(lowestindex < numvalues)
    {
        if(currentprime - startvalue < numvalues) lowestindex++;
        printf("sieve prime %d on %d.\n", currentprime, s);

        int offset = (startvalue/currentprime)*currentprime - startvalue;
        if (offset < 0)
            offset += currentprime;
        for(int i = offset; i < numvalues; i += currentprime)
        {
            values[i] = true;
        }
        long lowest = n + 1;
        for(int i = lowestindex; i < numvalues; i++)
        {
            if(values[i] == true) continue;
            lowest = startvalue + lowestindex;
            break;
        }
        printf("Lowest index is %d.\n", lowestindex);
        printf("Lowest is %d.\n", lowest);
        
        for (long t=0; t<p; t++)
        {
            printf("wrote %d to %d.\n", lowest, t);
            
            bsp_put(t,&lowest,Lowest,s*sizeof(long),sizeof(long));
        }
        bsp_sync();

        currentprime = lowest;
        for(long t = 0; t < p; t++)
        {
            if(Lowest[t] < currentprime)
            {
                currentprime = Lowest[t];
            }
            if(s == 1)
                printf("Lowest %d is %d.\n", t, Lowest[t]);
        }
        bsp_sync();

    }
    bsp_pop_reg(Lowest);
    free(Lowest);
    bsp_sync();


    bsp_end();
}

int main(int argc, char **argv){

    bsp_init(sieve, argc, argv);

    P=bsp_nprocs();
    if (argc>1) P=atol(argv[1]);
    if (P > bsp_nprocs()){
        printf("Sorry, only %u processors available.\n",
                bsp_nprocs());
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    N=1000;
    if (argc>2) N=atol(argv[2]);
    Num = (N+P)/P;

    /* SPMD part */
    sieve();

    /* Sequential part */
    exit(EXIT_SUCCESS);

} /* end main */