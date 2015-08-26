INCLUDE = /usr/local/include
#INCLUDE = ..
DBUGFLAGS = -g -traceback -check # debug version
FCFLAGS = -O3 -I$(INCLUDE) # release version
#LIBS = -lmpich -lpthread -lmpl
#FC = ifort
FC = /usr/local/bin/mpif90
# Linker flags (gfortran on OSX)
#LDFLAGS = -static-libgfortran -static-libgcc

TARGET = test_main

MODULES = inputs.o profile.o genome.o random_pkg.o init.o unittest.o sort.o\
			 mutation.o selection.o migration.o

OBJECTS = $(MODULES) $(TARGET).o fileio.o

all: $(OBJECTS) test_main.o
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJECTS) 

selection: $(MODULES) test_main.o
	$(FC) $(FCFLAGS) -o test_main $(OBJECTS)

f2py:
	f2py -c -m random_pkg random_pkg.f90
random:
	#f2py -c --fcompiler=intel --f90exec=/usr/bin/ifort -m random_pkg random_pkg.f90
	#f2py -c --f90exec=/usr/bin/ifort --f77exec=/usr/bin/ifort -m random_pkg random_pkg.f90

clean:
	rm *.o *.mod

test_main.o:	test_main.f90
	$(FC) $(FCFLAGS) -c test_main.f90
#test_selection.o:	test_selection.f90
#	$(FC) $(FCFLAGS) -c test_selection.f90
unittest.o:	unittest.f90
	$(FC) $(FCFLAGS) -c unittest.f90
inputs.o:	inputs.f90
	$(FC) $(FCFLAGS) -c inputs.f90
selection.o:    selection.f90 common.h genome.f90
	$(FC) $(FCFLAGS) -c selection.f90
random_pkg.o:  random_pkg.f90
	$(FC) $(FCFLAGS) -c random_pkg.f90
sort.o:     sort.f90
	$(FC) $(FCFLAGS) -c sort.f90
genome.o:       genome.f90
	$(FC) $(FCFLAGS) -c genome.f90
init.o:     init.f90 common.h
	$(FC) $(FCFLAGS) -c init.f90
profile.o:  profile.f90
	$(FC) $(FCFLAGS) -c profile.f90
mutation.o: mutation.f90
	$(FC) $(FCFLAGS) -c mutation.f90
migration.o:   migration.f90 common.h
	$(FC) $(FCFLAGS) -c migration.f90
fileio.o:   fileio.f90 common.h
	$(FC) $(FCFLAGS) -c fileio.f90
