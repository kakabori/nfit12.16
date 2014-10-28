_OBJS =	boxchain.o dataset.o funclmdif.o funsupport.o globalVariables.o \
        modelcalculator.o nrutil.o Para.o toad.o toadcmd.o toadmisc.o \
        tvds.o tvDSfit.o tvImg.o tvLinfitDriver.o utable.o \
        fileTools.o alglib_interpolation.o tcl_utility.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

CC = g++
ODIR = obj
SDIR = src
IDIR = inc

TCL = tcl8.4
TK = tk8.4
BLT = BLT

CFLAGS = -Wall -O3 -ffast-math -fPIC -g
CLINK = -lm -l$(TCL) -l$(TK) -ltiff -l$(BLT) -lpthread -lmydll \
        -lgsl -lgslcblas -Wl,-rpath,. -lalglib -linterp2d -lpthread
CPATH = -I/usr/include/$(TCL) -I/usr/local/include -I./$(IDIR) -I. -I../ \
        -L/usr/local/lib -L.

all: libtoad toad

toad: libtoad
	$(CC) $(ODIR)/toad.o $(CPATH) -ltoad_threaded $(CLINK) -o toad 

libtoad: $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) $(CPATH) $(CLINK) -shared -o libtoad_threaded.so 

$(ODIR)/%.o: $(SDIR)/%.cpp $(IDIR)/%.h
	$(CC) -c $(CPATH) $(CFLAGS) -o $@ $< 

clean:
	rm -f $(ODIR)/*.o libtoad_threaded.so


	

