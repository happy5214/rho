CC = gcc
LIBS = -lgmp
LOOP_COUNT = 1
MAX_ITERATIONS = 1000
FLAGS = -std=gnu99 -O2 -DLOOP_COUNT=$(LOOP_COUNT) -DMAX_ITERATIONS=$(MAX_ITERATIONS)

HEADERS = carg_parser.h factor.h rho.h types.h
objs = carg_parser.o rho.o io.o factor_common.o

floyd_objs = floyd.o $(objs)

.PHONY: all

all: floyd

floyd: $(floyd_objs)
	$(CC) -o $@ $(floyd_objs) $(LIBS)

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(FLAGS)

clean: 
	rm -f floyd *.o
