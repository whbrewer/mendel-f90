INCLUDE ?= /usr/local/include
SRC_DIR = ../src
MODDIR = ../src
PROJECT_INCLUDE = ../src
#INCLUDE = ..
LEGACYFLAGS ?= -fallow-argument-mismatch -std=legacy
DBUGFLAGS = -g -traceback -check $(LEGACYFLAGS) # debug version
FCFLAGS = -O3 -I$(PROJECT_INCLUDE) -I$(INCLUDE) -I$(MODDIR) -J$(MODDIR) $(LEGACYFLAGS) # release version
#LIBS = -lmpich -lpthread -lmpl
FC ?= mpif90
# Linker flags (gfortran on OSX)
#LDFLAGS = -static-libgfortran -static-libgcc

TARGET = test_main

MODULES = mpi_helpers.o inputs.o profile.o genome.o random_pkg.o polygenic.o init.o unittest.o sort.o\
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
	rm -f *.o *.mod

test_main.o:	test_main.f90
	$(FC) $(FCFLAGS) -c test_main.f90
#test_selection.o:	test_selection.f90
#	$(FC) $(FCFLAGS) -c test_selection.f90
unittest.o:	unittest.f90
	$(FC) $(FCFLAGS) -c unittest.f90
inputs.o:	$(SRC_DIR)/inputs.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/inputs.f90
selection.o:    $(SRC_DIR)/selection.f90 $(PROJECT_INCLUDE)/common.h $(SRC_DIR)/genome.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/selection.f90
random_pkg.o:  $(SRC_DIR)/random_pkg.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/random_pkg.f90
sort.o:     $(SRC_DIR)/sort.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/sort.f90
genome.o:       $(SRC_DIR)/genome.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/genome.f90
init.o:     $(SRC_DIR)/init.f90 $(PROJECT_INCLUDE)/common.h mpi_helpers.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/init.f90
mpi_helpers.o:	$(SRC_DIR)/mpi_helpers.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/mpi_helpers.f90
polygenic.o:	$(SRC_DIR)/polygenic.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/polygenic.f90
profile.o:  $(SRC_DIR)/profile.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/profile.f90
mutation.o: $(SRC_DIR)/mutation.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/mutation.f90
migration.o:   $(SRC_DIR)/migration.f90 $(PROJECT_INCLUDE)/common.h mpi_helpers.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/migration.f90
fileio.o:   $(SRC_DIR)/fileio.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/fileio.f90
