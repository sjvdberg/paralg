#include "bspedupack.h"

/*  This program computes the sum of the first n squares,
        sum = 1*1 + 2*2 + ... + n*n
    by computing the inner product of x=(1,2,...,n)^T and
    itself, for n>=0.
    The output should equal n*(n+1)*(2n+1)/6.
    The distribution of x is cyclic.
*/

static long P; // number of processors requested
static long N; // global vector length

/* this is the variant from the book: everyone sends their local
   result to everyone else, and the final sum is done redundantly.
 */
double bspip(long n, double *x, double *y)
{
    /* Compute inner product of vectors x and y
       of length n>=0 */

    long p= bsp_nprocs(); // p = number of processors obtained
    long s= bsp_pid();    // s = processor number

    double *Inprod= vecallocd(p);
    bsp_push_reg(Inprod,p*sizeof(double));
    bsp_sync();

    double inprod= 0.0;
    for (long i=0; i<nloc(p,s,n); i++)
        inprod += x[i]*y[i];

    for (long t=0; t<p; t++)
        bsp_put(t,&inprod,Inprod,s*sizeof(double),sizeof(double));

    bsp_sync();

    double alpha= 0.0;
    for (long t=0; t<p; t++)
        alpha += Inprod[t];

    bsp_pop_reg(Inprod);
    vecfreed(Inprod);

    return alpha;

} /* end bspip (allgather variant from the book) */

/* can you implement your own variant? Many options exist.

   Challenge: try to implement the "Butterfly reduction" algorithm:

        P0  P1  P2  P3
          \/  \/  \/
          /\  /\  /\
        P0  P1  P2  P3
          \____/
          /    \        (and similar for P1<->P3)
        P0  P1  P2  P3
*/
double my_bspip(long n, double *x, double *y)
{
    /* Compute inner product of vectors x and y
       of length n>=0 */

    long p= bsp_nprocs(); // p = number of processors obtained
    long s= bsp_pid();    // s = processor number

    double inprod= 0.0;
    for (long i=0; i<nloc(p,s,n); i++)
        inprod += x[i]*y[i];

    double alpha=-1.0;

    /* ... YOUR CODE GOES HERE ... */

    return alpha;
} /* end my_bspip */


void bspinprod()
{
    bsp_begin(P);
    long p= bsp_nprocs();
    long n = N;
    long s= bsp_pid();
    if (s==0){
        if(n<0)
            bsp_abort("Error in input: n is negative");
    }

    long nl= nloc(p,s,n);
    double *x= vecallocd(nl);
    for (long i=0; i<nl; i++){
        long iglob= i*p+s;
        x[i]= iglob+1;
    }

    bsp_sync();
    double alpha;
    for (int type=0; type<=1; type++)
    {
      double time0= bsp_time();

      for (int run=0; run<100; run++)
      {
        if (type==0)
        {
          alpha = bspip(n,x,x); // original implementation in BSPedupack
        }
        else if (type==1)
        {
          alpha = my_bspip(n,x,x);
        }
        else
        {
          bsp_abort("'type' got an invalid value, must be 0 (allgather), 1 (gather/scatter) or 2 (butterfly)");
        }
      }
      bsp_sync();
      double time1= bsp_time();

      printf("Proc %ld: sum of squares up to %ld*%ld is %.lf\n",
            s,n,n,alpha); fflush(stdout);
      if (s==0){
          printf("np=%d, average time for variant BSP-%d: %.6lf seconds.\n", p, type, (time1-time0)/100.0);
      }
      bsp_sync();
    }

    /* Compute exact output (number can become large) */
    double sum = (n*(n+1.0)*(2.0*n+1)) / 6.0;
    printf("n(n+1)(2n+1)/6 = %lf\n", sum); fflush(stdout);

    vecfreed(x);
    bsp_end();

} /* end bspinprod */

// main program to benchmark the inner product s=x^T*y in three different implementations.

// Usage: ./bsp_inprod.x <p> <N>
// where: p is the number of processes/threads to be used,
//        N is the global vector length.
int main(int argc, char **argv){

    bsp_init(bspinprod, argc, argv);

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
