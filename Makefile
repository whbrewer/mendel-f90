INSTALL_DIR = /usr/local/bin
SRC_DIR = src
PROJECT_INCLUDE = include

GIT_VERSION := $(shell git describe --abbrev=7 --dirty --always --tags)

#FC = /opt/intel/fc/10.0.026/bin/ifort -vec-report0
#FC = /opt/intel/bin/ifort
#FC = /opt/pgi/linux86-64/8.0-4/bin/pgf90 # c101
#FC = /usr/local/bin/mpif90
FC = gfortran

# Following are needed for building parallel version
# Comment out if compiling with mpif90
# or if making serial version.
# MPICH
#LIBS = -lmpich -lpthread -lmpl
#LIBS = -L/usr/local/lib -lmpich -lpthread -lmpl
# Open MPI
#LIBS = -L/usr/lib64/mpi/gcc/openmpi/lib64 -lopen-rte -lmpi

# when using Open MPI
#INCLUDE = /usr/lib64/mpi/gcc/openmpi/include
# when using MPICH
INCLUDE = /usr/local/lib

INCLUDE = /usr/local/include
# Compiler flags
DBUGFLAGS = -g -traceback -check # debug version
#FCFLAGS = -traceback -O3 -I$(INCLUDE) # release version ifort
FCFLAGS = -O3 -I$(PROJECT_INCLUDE) -I$(INCLUDE) # release version gfortran
# note use flag -fpe:0 to handle floating point exceptions

# Linker flags (gfortran on OSX)
#LDFLAGS = -static-libgfortran -static-libgcc

SERIALFN = mendel_serial

# executable name
TARGET = mendel

MODULES = sort.o random_pkg.o inputs.o genome.o profile.o polygenic.o \
          init.o selection.o

OTHERS = $(MODULES) mutation.o mating.o fileio.o

POBJECTS = $(OTHERS) diagnostics.o mendel.o migration.o

SERIAL_OTHERS = $(filter-out init.o,$(OTHERS))
SOBJECTS = $(SERIAL_OTHERS) $(SERIALFN).o init_serial.o serial_stubs.o

TOBJECTS = $(OTHERS) diagnostics.o test.o migration.o

##########################################
# build rules
##########################################

all: release

debug: FCFLAGS = $(DBUGFLAGS)
debug: $(POBJECTS)
	$(FC) $(DBUGFLAGS) $(LDFLAGS) -o $(TARGET) $(POBJECTS) $(LIBS)

test: $(TOBJECTS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o test $(TOBJECTS) $(LIBS)

pre-build:
	sed -i.bak 's/VERSION.*VERSION/VERSION >>> $(GIT_VERSION) <<< VERSION/' $(SRC_DIR)/init.f90

release: pre-build $(POBJECTS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $(TARGET) $(POBJECTS) $(LIBS)

# in order to build the serial version, first
# run "make preserial", then "make serial"
parallel: release

serial: $(SOBJECTS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $(SERIALFN) $(SOBJECTS) $(LIBS)

install:
	install $(TARGET) $(INSTALL_DIR)

uninstall:
	echo "removing file $(INSTALL_DIR)/$(TARGET)"
	rm $(INSTALL_DIR)/$(TARGET)

dist:
	scripts/package.sh

preserial:
	cp $(SRC_DIR)/mendel.f90 $(SRC_DIR)/$(SERIALFN).f90
	cat $(SRC_DIR)/diagnostics.f90 >> $(SRC_DIR)/$(SERIALFN).f90
	cp $(SRC_DIR)/init.f90 $(SRC_DIR)/init_serial.f90
	sed -i.bak '/START_MPI/,/END_MPI/d' $(SRC_DIR)/init_serial.f90
	sed -i.bak '/START_MPI/,/END_MPI/d' $(SRC_DIR)/$(SERIALFN).f90

cln:
	\rm -f mendel.o migration.o mendel

clean:
	\rm -f *.o *.mod $(TARGET) test0* $(SRC_DIR)/*_serial.f90 *.f90-e a.out\
	       success mendel_serial

###########################################
# dependencies
###########################################

sort.o:		$(SRC_DIR)/sort.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/sort.f90

random_pkg.o:	$(SRC_DIR)/random_pkg.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/random_pkg.f90

mendel.o:       $(SRC_DIR)/mendel.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/mendel.f90

init.o:		$(SRC_DIR)/init.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/init.f90

init_serial.o:		$(SRC_DIR)/init_serial.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/init_serial.f90

serial_stubs.o:	serial_stubs.f90
	$(FC) $(FCFLAGS) -c serial_stubs.f90

diagnostics.o:  $(SRC_DIR)/diagnostics.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/diagnostics.f90

fileio.o:	$(SRC_DIR)/fileio.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/fileio.f90

selection.o:    $(SRC_DIR)/selection.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/selection.f90

mating.o:	$(SRC_DIR)/mating.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/mating.f90

mutation.o:	$(SRC_DIR)/mutation.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/mutation.f90

migration.o:	$(SRC_DIR)/migration.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/migration.f90

mendel_serial.o:	$(SRC_DIR)/mendel_serial.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/mendel_serial.f90

migration_serial.o:	$(SRC_DIR)/migration_serial.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/migration_serial.f90

polygenic.o:	$(SRC_DIR)/polygenic.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/polygenic.f90

test.o:		$(SRC_DIR)/test.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/test.f90

profile.o:	$(SRC_DIR)/profile.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/profile.f90

inputs.o:       $(SRC_DIR)/inputs.f90
	        $(FC) $(FCFLAGS) -c $(SRC_DIR)/inputs.f90

genome.o:       $(SRC_DIR)/genome.f90
	        $(FC) $(FCFLAGS) -c $(SRC_DIR)/genome.f90
