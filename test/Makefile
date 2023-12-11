
# for debugging/development, sue:
#CFLAGS=-O0 -g -fsanitize=bounds -march=native
# for benchmarking/performance evaluation, use:
CFLAGS=-O3 -march=native

default: bsp

bsp: test.c


test.x: test.c
	bspcc ${CFLAGS} -o test.x 

clean:
	-rm *.o *.x
