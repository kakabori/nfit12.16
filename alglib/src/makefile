CXX = g++
CXXFLAGS = 
CXXPATH = -I.
CXXLINK = 
DEPS = \
	alglibinternal.h \
	alglibmisc.h \
	ap.h \
	dataanalysis.h \
	diffequations.h \
	fasttransforms.h \
	integration.h \
	interpolation.h \
	linalg.h \
	optimization.h \
	solvers.h \
	specialfunctions.h \
	statistics.h \
	stdafx.h
OBJ = \
	alglibinternal.o \
	alglibmisc.o \
	ap.o \
	dataanalysis.o \
	diffequations.o \
	fasttransforms.o \
	integration.o \
	interpolation.o \
	linalg.o \
	optimization.o \
	solvers.o \
	specialfunctions.o \
	statistics.o

all: $(DEPS) libalglib

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXPATH)

libalglib: $(OBJ)
	$(CXX) $(OBJ) $(CXXFLAGS) $(CXXPATH) $(CXXLINK) -shared -o libalglib.so 

clean:
	rm -f *.o
