_OBJS =	alglib_interpolation.o test_Alglib_CubicSpline2D.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

CC = g++
ODIR = ../obj
SDIR = ../src
IDIR = ../inc

#CFLAGS = -Wall -O0 -ffast-math -fPIC -lpthread -g
CFLAGS = -Wall -O3 -ffast-math -fPIC -lpthread
CLINK = -lm -lpthread -lalglib -lgsl -lgslcblas -Wl,-rpath,.
CPATH = -I$(IDIR) -I../../ -L../ -I../

all: test_Alglib_CubicSpline2D_class

test_Alglib_CubicSpline2D_class: $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) $(CPATH) $(CLINK) -o test_Alglib_CubicSpline2D_class

$(ODIR)/%.o: $(SDIR)/%.cpp $(IDIR)/%.h
	$(CC) -c $(CPATH) $(CFLAGS) -o $@ $< 

clean:
	rm -f $(ODIR)/test_Alglib_CubicSpline2D.o \
	      $(ODIR)/alglib_interpolation.o


	

