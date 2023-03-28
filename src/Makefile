## Adapted from https://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

IDIR = ../include
CXX = g++
CFLAGS = -std=c++11 -Wall -g -I$(IDIR)

ODIR = obj

LIBS=-lm

_DEPS = Anfangsverteilung.h MersenneTwister.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = ecDNA_evolution.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

ecDNA_evolution: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LIBS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) 

.PHONY : clean

clean:
	rm -f $(ODIR)/*.o ecDNA_evolution
