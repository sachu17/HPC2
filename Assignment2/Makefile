#  Sample makefile for the Henon basin boundary problem

SOURCES = precision.f90 henondim.f90 main.f90
OBJECTS  = precision.o henondim.o main.o

#
# Gfortran on the Lab computers
#
F95 = gfortran
LIBS = 
PARALLEL = 
DEBUG = -g 
OPT = -O2 -ftree-vectorize
OTHER = -c -Wall 

#-----------------------------------------------------------------------
#  Choose one or the other
FFLAGS = $(OPT) $(OTHER)
#FFLAGS = $(DEBUG) $(OTHER)
#-----------------------------------------------------------------------

%.o: %.f90
	$(F95)  $(FFLAGS) $<
henondim: $(OBJECTS)
	$(F95) $(PARALLEL) $(OBJECTS) $(LIBS) -o henondim

precision.o: 
henondim.o: precision.o
main.o: henondim.o

clean:
	rm -f $(OBJECTS) a.out core *.mod
