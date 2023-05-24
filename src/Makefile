## Adapted from https://www.cs.colby.edu/maxwell/courses/tutorials/maketutor/

IDIR = ../include
CXX = gcc
CFLAGS = -std=c++11 -Wall -g -I$(IDIR)

ODIR = obj

LIBS=-lm -lboost_filesystem -lboost_system -lstdc++ 

_DEPS = Anfangsverteilung.h MersenneTwister.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = ecDNA_evolution_many_types.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

ecDNA_evolution_many_types.exe: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS) $(LIBS)

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS) 

.PHONY : clean

clean:
	rm -f $(ODIR)/*.o *.txt ecDNA_evolution_many_types.exe
