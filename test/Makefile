
# for debugging/development, sue:
#CFLAGS=-O0 -g -fsanitize=bounds -march=native
# for benchmarking/performance evaluation, use:
CFLAGS=-O3 -march=native

default: bsp

bsp: bsp_inprod.x


compiled.x: bsp_inprod.c
	bspcc ${CFLAGS} -o compiled.x 

clean:
	-rm *.o *.x
