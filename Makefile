#*******************************************************************************
# This file is part of GadgetSMP.
# Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
#******************************************************************************/



CC=gcc
FC=gfortran

OPTS=-DTIMER -DCOUNTER -DCPU_TPS=2.67e9 -DHAVE_ZLIB -DHAVE_SETAFFINITY -D_GNU_SOURCE
CFLAGS=-O3 -g -std=gnu99 -Wall -Werror -march=native -mtune=native -ffast-math -fomit-frame-pointer -malign-double -fstrict-aliasing -fopenmp
# CFLAGS=-O0 -g -std=gnu99 -Wall -Werror -fopenmp
LDFLAGS=-std=gnu99 -lm -lpthread -fopenmp -lz

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
