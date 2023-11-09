
# for debugging/development, sue:
#CFLAGS=-O0 -g -fsanitize=bounds -march=native
# for benchmarking/performance evaluation, use:
CFLAGS=-O3 -march=native

default: bsp

bsp: parallel.x


bsp_inprod.x: parallel.c bspedupack.c
	bspcc ${CFLAGS} -o parallel.x $< bspedupack.c -lm


clean:
	-rm *.o *.x
