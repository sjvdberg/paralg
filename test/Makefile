
# for debugging/development, sue:
#CFLAGS=-O0 -g -fsanitize=bounds -march=native
# for benchmarking/performance evaluation, use:
CFLAGS=-O3 -march=native

default: mpi

mpi: test.x
bsp: test.x


test.x: test.c 
	~/MultiCoreBSP-for-C/tools/bspcc ${CFLAGS} -o bsp_inprod.x $< bspedupack.c -lm

test.x: test.c
	mpicc ${CFLAGS} -o test.x $< -lm

clean:
	-rm *.o *.x