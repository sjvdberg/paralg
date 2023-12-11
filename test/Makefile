
# for debugging/development, sue:
#CFLAGS=-O0 -g -fsanitize=bounds -march=native
# for benchmarking/performance evaluation, use:
CFLAGS=-O3 -march=native

default: bsp

bsp: code.c


compiled.x: code.c
	bspcc ${CFLAGS} -o compiled.x 

clean:
	-rm *.o *.x
