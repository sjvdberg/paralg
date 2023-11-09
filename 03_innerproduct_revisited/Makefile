
# for debugging/development, sue:
#CFLAGS=-O0 -g -fsanitize=bounds -march=native
# for benchmarking/performance evaluation, use:
CFLAGS=-O3 -march=native

default: mpi

mpi: mpi_inprod.x
bsp: bsp_inprod.x


bsp_inprod.x: bsp_inprod.c bspedupack.c
	bspcc ${CFLAGS} -o bsp_inprod.x $< bspedupack.c -lm

mpi_inprod.x: mpi_inprod.c
	mpicc ${CFLAGS} -o mpi_inprod.x $< -lm

clean:
	-rm *.o *.x
