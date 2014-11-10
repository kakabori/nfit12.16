_OBJS =	funsupport.o globalVariables.o BareStructureFactor.o Para.o \
        alglib_interpolation.o tvds.o utable.o fileTools.o \
        create_structure_factor.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

CC = g++
ODIR = ../obj
SDIR = ../src
IDIR = ../inc

TCL = tcl8.4
TK = tk8.4
BLT = BLT

#CFLAGS = -Wall -O0 -ffast-math -fPIC -lpthread -g
CFLAGS = -Wall -O3 -ffast-math -fPIC -std=c++0x
CLINK = -lm -l$(TCL) -l$(TK) -l$(BLT) -lpthread -lalglib\
        -lgsl -lgslcblas -Wl,-rpath,.
CPATH = -I/usr/include/$(TCL) -I$(IDIR) -I../../ -L../ -I../

all: create_structure_factor

create_structure_factor: $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) $(CPATH) $(CLINK) -o create_structure_factor

$(ODIR)/%.o: $(SDIR)/%.cpp $(IDIR)/%.h
	$(CC) -c $(CPATH) $(CFLAGS) -o $@ $< 

clean:
	rm -f $(ODIR)/funsupport.o $(ODIR)/globalVariables.o \
	      $(ODIR)/BareStructureFactor.o $(ODIR)/Para.o $(ODIR)/tvds.o \
	      $(ODIR)/utable.o $(ODIR)/fileTools.o \
	      $(ODIR)/create_structure_factor.o \
	      $(ODIR)/alglib_interpolation.o

