
CC=gcc
FC=gfortran

OPTS=-DTIMER -DCOUNTER -DCPU_TPS=2.67e9
CFLAGS=-O3 -g -std=gnu99 -Wall -Werror -march=native -mtune=native -ffast-math -fomit-frame-pointer -malign-double -fstrict-aliasing -fopenmp
# CFLAGS=-O0 -g -std=gnu99 -Wall -Werror -fopenmp
LDFLAGS=-lm -lpthread -fopenmp

FFLAGS=$(CFLAGS)

OBJS=space.o runner.o test.o

default: all

depfile: $(wildcard *.c) $(wildcard *.h)
	$(CC) -MM *.c *.h > depfile
    
include depfile

%.o: Makefile

%.o: %.c Makefile
	$(CC) -c $(CFLAGS) $(OPTS) $< -o $@
# 	$(CC) -c $(CFLAGS) $(OPTS) $< -dA -S

test: $(OBJS)
	gcc $^ $(LDFLAGS) -o test

all: test 

clean:
	rm -f $(OBJS) test
