SOURCES = relaxation.f90 main.f90
OBJECTS  = relaxation.o main.o 

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

relaxation: $(OBJECTS)
	$(F95) $(PARALLEL) $(OBJECTS) $(LIBS) -o relaxation


relaxation.o: 

main.o: relaxation.o


clean:
	rm -f $(OBJECTS) a.out core *.mod
