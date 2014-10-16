C = gcc
CPATH = -I.
DEPS = \
	clencurt.h \
	converged.h \
	cubature.h \
	vwrapper.h
OBJ = \
	clencurt_gen.o \
	hcubature.o \
	pcubature.o 

all: $(DEPS) libcubature

%.o: %.c $(DEPS)
	$(C) -c -o $@ $< $(CPATH)

libcubature: $(OBJ)
	$(C) $(OBJ) $(CPATH) -shared -o libcubature.so 

clean:
	rm -f *.o
