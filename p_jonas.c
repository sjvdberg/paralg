#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef unsigned long long int int_type;
static const MPI_Datatype mpi_int_type=MPI_UNSIGNED_LONG_LONG;

static inline int_type min(int_type a, int_type b) {return a<b? a: b;}
static inline int_type max(int_type a, int_type b) {return a>b? a: b;}

// this is a sequential function that takes a list of num_primes
// prime numbers 'primes' and filters the range [smin, smax] to extend the
// list. Returns the number of primes found in the interval [smin,smax],
// and the list of new prime numbers in new_primes. On entry, new_primes must
// have space allocated for at least (smax-smin)/2 + 1 numbers.
int_type sieve(int_type num_primes_in, int_type* primes, int_type* new_primes, int_type smin, int_type smax)
{
  int_type chunk_size=max(100000, smax/100);

  int_type *test = (int_type*)malloc(chunk_size*sizeof(int_type));

  int_type new_num_primes=0;

  for (int_type ic=smin; ic<=smax; ic+=chunk_size)
  {
    int_type imax=min(smax, ic+chunk_size-1);
    for (int_type i=ic; i<=imax; i++) test[i-ic] = i;
    for (int_type j=0; j<num_primes_in; j++)
    {
      int_type p=primes[j];
      if (imax<p*p) break;
      int_type imin=ic-ic%p; if (imin<ic) imin+=p;
      for (int_type i=imin; i<=imax; i+=p)
      {
        test[i-ic]=0;
      }
    }
    for (int_type s=ic; s<=imax; s++)
    {
      if (test[s-ic])
      {
        new_primes[new_num_primes++] = test[s-ic];
      }
    }
  }
/*
  if (new_num_primes==0) fprintf(stdout, " -> no new primes found\n");
  if (new_num_primes==1) fprintf(stdout, " -> new prime: [%lld]\n", new_primes[0]);
  if (new_num_primes==2) fprintf(stdout, " -> %lld new primes: [%lld, %lld]\n", new_num_primes, new_primes[0], new_primes[new_num_primes-1]);
  if (new_num_primes>2) fprintf(stdout, " -> %lld new primes: [%lld, ..., %lld]\n", new_num_primes, new_primes[0], new_primes[new_num_primes-1]);
*/
  free(test);
  return new_num_primes;
}

// Estimate the number of primes <=N
// See octave's "help primes":
//
// The distance from one prime to the next is, on average,
// proportional to the logarithm of the prime.  Integrating, one finds
// that there are about k primes less than k*log (5*k).
int_type nprimes(int_type N)
{
  if (N==1) return 0;
  if (N==2) return 1;
  double k = (double)N/2.0;
  double kprev=0;
  // Newton iteration
  for (int i=0; i<5; i++)
  {
    kprev=k;
    double logk=log(k);
    double f = N - k*logk;
    double fp= 1 - (logk + 1);
    k=kprev-f/fp;
  }
  return (int_type)k;
}

int_type calc_primes(int_type N, int_type* primes)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, nproc;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nproc);

  int *new_num_primes = (int*)malloc(nproc*sizeof(int));
  int *new_prime_disps = (int*)malloc(nproc*sizeof(int));

  primes[0]=2;
  int_type num_primes=1;
  int_type smin, smax=min(N,3);

  int max_my_nprimes = 0;
  int_type* my_primes = NULL;

  while (1)
  {
    smin=primes[num_primes-1]+1;
    smax=min(smax, N);
    // sequential version:
    //num_primes += sieve(num_primes, primes, primes+num_primes, smin, smax);
    int_type sdif = (smax-smin+1)/nproc;
    int_type my_smin = smin+rank*sdif;
    int_type my_smax = my_smin+sdif;
    if (rank==nproc-1) my_smax=smax;

    fprintf(stdout, "P%d, smin=%lld, smax=%lld\n",rank,my_smin, my_smax);

    int_type my_nprimes = nprimes(my_smax)-nprimes(my_smin)+1;
    if (my_nprimes>max_my_nprimes)
    {
      if (primes) free(my_primes);
      max_my_nprimes=my_nprimes;
      my_primes = (int_type*)malloc(my_nprimes*sizeof(int_type));
    }

    int_type my_new_num_primes = sieve(num_primes, primes, my_primes, my_smin, my_smax);
    //Allgather the number of new primes found by each process in their subinterval
    MPI_Allgather(&my_new_num_primes, 1, mpi_int_type, new_num_primes, 1, mpi_int_type, comm);
    new_prime_disps[0]=0;
    for (int i=0; i<nproc-1; i++) new_prime_disps[i+1]=new_prime_disps[i]+new_num_primes[i];
    // Allgather the new prime numbers from this interval
    MPI_Allgatherv(my_primes, my_new_num_primes, mpi_int_type, primes+num_primes, new_num_primes, new_prime_disps, mpi_int_type, comm);
    for (int i=0; i<nproc; i++) num_primes+=new_num_primes[i];
    if (smax==N) break;
    smax=smax*smax;
  }
  if (my_primes) free(my_primes);
  free(new_num_primes);
  free(new_prime_disps);
  return num_primes;
}

// runs a number of tests with known outcome and returns -f, where f is the number of failed tests.
int run_tests(int rank)
{
  const int ntests  = 6;
  int npassed       = 0;
  const int_type Ns[6]                = {2, 3, 13, 100, 1000, 1000000000};
  const int_type num_primes[6]   = {1, 2, 6,   25,  168, 50847534};
  const int_type max_prime[6]    = {2, 3, 13,  97,  997, 999999937};

  int_type *primes = (int_type*)malloc(Ns[ntests-1]*sizeof(int_type));

  for (int i=0; i<ntests; i++)
  {
    int_type N=Ns[i];
    int_type np = calc_primes(N, primes);
    if (np<=nprimes(N) && np==num_primes[i] && primes[np-1]==max_prime[i])
    {
      npassed++;
    }
    else if (rank==0)
    {
      fprintf(stdout, "TEST FAILED: N=%lld\n", N);
      fprintf(stdout, "             expect at most %lld primes in interval, found %lld\n", nprimes(N), np);
      fprintf(stdout, "             found %lld primes, expected %lld\n", np, num_primes[i]);
      fprintf(stdout, "             max prime is %lld, expected %lld\n", primes[np-1], max_prime[i]);
    }
  }
  return npassed-ntests;
}

// MPI program to compute all prime numbers smaller or equal to N,
// where N is read from the command line: e.g.,
//
// mpirun -np 4 ./sieve.x 1000
//
// computes [2, 3, 5, ..., 997] (stored distributed among all processes),
// and prints out the largest prime number found, and the time it took to
// compute all these prime numbers.
//
// The memory requirement can be estimated based on the following observation
// (see octave's "help primes"):
//
// The distance from one prime to the next is, on average,
// proportional to the logarithm of the prime.  Integrating, one finds
// that there are about k primes less than k*log (5*k).
//
// Consequently, we expect about


// If no argument is given, we run a couple of tests and print SUCCESS or FAILURE at the end.
//
int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int_type N=0;
  if (argc>1) N = atoll(argv[1]);

  int rank=0, nproc=1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (N==0)
  {
    if (rank==0)
    {
      fprintf(stdout,"no argument N given -- run tests.\n");
    }
    int result=run_tests(rank);
    if (rank==0)
    {
      if (result==0)
      {
        fprintf(stdout, "ALL TESTS PASSED\n");
      }
      else
      {
        fprintf(stdout, "%d TESTS FAILED\n", -result);
      }
    }
    MPI_Finalize();
    return result;
  }

  int_type k = nprimes(N);
  fprintf(stdout, "Expecting k=%lld primes in [0, %lld]\n", k, N);

  int_type *primes = (int_type*)malloc(k*sizeof(int_type));

  if (primes==NULL)
  {
    if (rank==0) fprintf(stderr, "Memory allocation failed - I need about %e MBytes per process for this calculation.\n", (double)k*sizeof(int_type)*8.0e-6);
    MPI_Abort(MPI_COMM_WORLD,-1);
  }

  double t0 = MPI_Wtime();

  int_type num_primes=calc_primes(N, primes);

  double t1 = MPI_Wtime();

  int_type max_prime = primes[num_primes-1];
  if (rank==0)
  {
    fprintf(stdout, "Found %lld prime numbers in the interval [2, %lld].\n", num_primes, N);
    fprintf(stdout, "Largest prime number found: %lld\n", primes[num_primes-1]);
    fprintf(stdout, "Runtime: %e seconds.\n", t1-t0);
    if (N<=100)
    {
      fprintf(stdout, "Primes found: \n");
      for (int i=0; i<num_primes; i++) fprintf(stdout, "%lld\n", primes[i]);
    }
  }
  MPI_Finalize();
  free(primes);
  return 0;
}
