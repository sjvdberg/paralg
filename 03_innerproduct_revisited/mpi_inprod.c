#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

/*  This program computes the sum of the first n squares, for n>=0,
        sum = 1*1 + 2*2 + ... + n*n
    by computing the inner product of x=(1,2,...,n)^T and itself.
    The output should equal n*(n+1)*(2n+1)/6.
    The distribution of x is cyclic.
*/

// useful macro to check the return value from MPI calls
#define CHECK_MPI(s) {\
  int iflag=(s);\
  if (iflag!=MPI_SUCCESS) {\
    fprintf(stderr, "MPI call '%s' returned value %d. Aborting -- check the corresponding man-page for the meaning of the error.\n",\
    (#s), iflag);\
    MPI_Abort(MPI_COMM_WORLD, iflag);\
  }\
}

double my_allreduce(double alpha, MPI_Comm comm, int variant)
{
  int np, rank;
  double result=0.0;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &np);

  if (variant==0)
  {
    // allgather + local sum, as implemented in the BSPedupack inprod driver
    double tmp[np];
    // note: we use non-blocking send and receives because otherwise the order
    // will matter and a deadlock may occur.
    MPI_Request requests[2*np];
    for (int p=0; p<np; p++)
    {
      // send local value to rank p and mark it with a tag corresponding to my rank.
      // note that this tag is arbitrary, but the receiving process has to use the same
      // tag when posting his Irecv-call.
      int send_tag=rank;
      CHECK_MPI(MPI_Isend(&alpha, 1, MPI_DOUBLE, p, send_tag, comm, &requests[p]));
      // issue the receive call to get the partial result from rank p.
      // Note that they will have used the tag p to send it.
      int recv_tag=p;
      CHECK_MPI(MPI_Irecv(tmp+p,1,MPI_DOUBLE,p,recv_tag,comm,&requests[np+p]));
    }
    // Now make sure that all communication has finished and compute the overall sum.
    // Note that no additional barrier is required!
    CHECK_MPI(MPI_Waitall(2*np, requests,MPI_STATUSES_IGNORE));
    for (int i=0; i<np; i++) result+=tmp[i];
  }
  else if (variant==1)
  {
    /* ... YOUR CODE GOES HERE ... */
    result = 0.0; /* should be the sum of all alpha's from the processes */
  }
  else
  {
    // default: MPI_Allreduce, should always be the fastest variant
    CHECK_MPI(MPI_Allreduce(&alpha,&result,1,MPI_DOUBLE,MPI_SUM,comm));
  }
  return result;
}

static int nloc(int p, int s, int n){
    /* Compute number of local components of processor s for vector
       of length n distributed cyclically over p processors. */

    return  (n+p-s-1)/p ;

} /* end nloc */

double mpiip(int p, int s, int n, double *x, double *y, int variant){
    /* Compute inner product of vectors x and y of length n>=0 */

    int nloc(int p, int s, int n);
    double alpha;
    int i;

    alpha = 0.0;
    int nl=nloc(p,s,n);
    for (i=0; i<nl; i++){
        alpha += x[i]*y[i];
    }
    alpha = my_allreduce(alpha, MPI_COMM_WORLD, variant);
    return alpha;

} /* end mpiip */


/* main program. Takes the vector length n and variant to use as command-line
   arguments, for example: ./mpi_inprod.x 100 0
   The variants are:
   0: manual implementation of "allgather", as in the book
   1: your own implementation
   any other: MPI_Allreduce
*/
int main(int argc, char **argv)
{
    double *x, alpha, time0, time1;
    int p, s;

    /* SPMD part */
    MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD,&p); /* p = number of processors */ 
    MPI_Comm_rank(MPI_COMM_WORLD,&s); /* s = processor number */ 

    int n=1000;
    if (argc>1) n=atoi(argv[1]);


    if(n<0) MPI_Abort(MPI_COMM_WORLD,-1);

    int nl= nloc(p,s,n);
    x= (double*)malloc(nl*sizeof(double));
    for (int i=0; i<nl; i++) x[i]= i*p+s+1;

    for (int variant=0; variant<=2; variant++)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      time0=MPI_Wtime();
      for (int run=0; run<100; run++)
      {
        alpha= mpiip(p,s,n,x,x,variant);
      }
      /* use barrier to get consistent timing across processes
         Note that there are no other barriers required in the program.
       */
      MPI_Barrier(MPI_COMM_WORLD); 
      time1=MPI_Wtime();
      printf("Processor %d, variant %d: sum of squares up to %d*%d is %.lf\n",
                s,variant, n,n,alpha); fflush(stdout);
      if (s==0)
      {
        printf("np=%d, average time variant MPI-%d: %.6le seconds\n",p, variant, (time1-time0)/100.0);
      }
    }
    free(x);
    MPI_Finalize();
    return 0;
} /* end main */
