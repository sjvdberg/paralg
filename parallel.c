#include <stdio.h>
#include <stdbool.h>
#include <math.h>

static long P; // number of processors requested
static long N; // upper bound on the primes.

void sieve()
{
    bsp_begin(P);
    long p = P;
    long s = bsp_pid();
    long n = N;

    int numvalues, startvalue;
    numvalues = (n + p) / p;
    startvalue = s * numvalues ;

    bool values[numvalues] = {false};
    int lowestindex = 0;
    int currentprime = 2;
    if(s == 0)
    {
        values[0] = true;
        values[1] = true;
        lowestindex = 2;
    }
    double *Lowest= malloc(MAX(n,1)*sizeof(double complex));
    bsp_push_reg(Lowest,p*sizeof(long));
    bsp_sync();

    while(lowestindex < numvalues)
    {
        if(s == 0)
            printtf("%d", currentprime)

        int offset = (startvalue/currentprime) - startvalue;
        if (offset < 0)
            offset += currentprime;
        if(start < startvalue)
            start += currentprime;
        for(int i = offset; i < numvalues; i += currentprime)
        {
            values[i] = true;
        }
        int lowest = n + 1;
        for(i = lowestindex; i < numvalues; i++)
        {
            if(values[i] == True) continue;
            lowest = startvalue + lowestindex
            break;
        }

        for (long t=0; t<p; t++)
            bsp_put(t,&lowest,Lowest,s*sizeof(long),sizeof(long));
        bsp.sync();

        currentprime = lowest
        for(long t = 0; t < p; t++)
            if(Lowest[t] < currentprime)
            {
                currentprime = Lowest[t];
            }
    }
    bsp.sync()
    bsp_pop_reg(Lowest);
    free(Lowest)

    bsp_end()
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

    /* SPMD part */
    bspinprod();

    /* Sequential part */
    exit(EXIT_SUCCESS);

} /* end main */