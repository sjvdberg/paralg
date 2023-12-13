
# for debugging/development, sue:
#CFLAGS=-O0 -g -fsanitize=bounds -march=native
# for benchmarking/performance evaluation, use:
CFLAGS=-O3 -march=native

default: mpi

mpi: parallel.x


parallel.x: parallel.c
	mpicc ${CFLAGS} -o parallel.x $< -lm

clean:
	-rm *.o *.x