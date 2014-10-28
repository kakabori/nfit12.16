CC = gcc

CFLAGS = -Wall -O3 -fopenmp
CLINK = -lm -lcubature -lstdc++ -lfftw3 -lfftw3l
CPATH = -I../cubature-1.0/ -L. -Wl,-rpath,. 

all: test_cubature

test_cubature: test_cubature.cpp
	$(CC) test_cubature.cpp $(CFLAGS) $(CPATH) $(CLINK) -o test_cubature

clean:
	rm -f test_cubature
