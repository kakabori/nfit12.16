_OBJS =	funsupport.o globalVariables.o modelcalculator.o Para.o \
        tvds.o utable.o fileTools.o test_ModelCalculator.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

CC = g++
ODIR = ../obj
SDIR = ../src
IDIR = ../inc

TCL = tcl8.4
TK = tk8.4
BLT = BLT

CFLAGS = -Wall -O0 -ffast-math -fPIC -lpthread -g
CLINK = -lm -l$(TCL) -l$(TK) -l$(BLT) -lpthread \
        -lgsl -lgslcblas -Wl,-rpath,.
CPATH = -I/usr/include/$(TCL) -I$(IDIR) 

all: test_ModelCalculator_class

test_ModelCalculator_class: $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) $(CPATH) $(CLINK) -o test_ModelCalculator_class

$(ODIR)/%.o: $(SDIR)/%.cpp $(IDIR)/%.h
	$(CC) -c $(CPATH) $(CFLAGS) -o $@ $< 

clean:
	rm -f $(ODIR)/funsupport.o $(ODIR)/globalVariables.o \
	      $(ODIR)/modelcalculator.o $(ODIR)/Para.o $(ODIR)/tvds.o \
	      $(ODIR)/utable.o $(ODIR)/fileTools.o \
	      $(ODIR)/test_ModelCalculator.o


	

