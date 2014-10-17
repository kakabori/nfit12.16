sources = spline2d_bicubic.cpp
CC = g++
CPATH = -I../src -L/usr/local/lib
CLINK = -lalglib

all: $(sources)
	$(CC) $^ $(CPATH) $(CLINK)

clean:
	rm -f *~
	rm -f *.o
	rm -f *#
	rm -f *.d
	rm -f *.d.*
	rm -f a.out
	rm -f *.dat
