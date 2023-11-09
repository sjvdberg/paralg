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
    printf("start sieve with %d processors.\n", p);


    while(lowestindex < numvalues)
    {
        printf("sieve prime %d\n", currentprime);

        int offset = (startvalue/currentprime) - startvalue;
        if (offset < 0)
            offset += currentprime;
        for(int i = offset; i < numvalues; i += currentprime)
        {
            values[i] = true;
        }
        int lowest = n + 1;
        for(int i = lowestindex; i < numvalues; i++)
        {
            if(values[i] == true) continue;
            lowest = startvalue + lowestindex;
            break;
        }

        for (long t=0; t<p; t++)
            bsp_put(t,&lowest,Lowest,s*sizeof(long),sizeof(long));
        bsp_sync();

        currentprime = lowest;
        for(long t = 0; t < p; t++)
            if(Lowest[t] < currentprime)
            {
                currentprime = Lowest[t];
            }

    }
    bsp_sync();
    bsp_pop_reg(Lowest);
    free(Lowest);

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