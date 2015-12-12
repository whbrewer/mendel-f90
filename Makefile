INSTALL_DIR = /usr/local/bin

#FC = /opt/intel/fc/10.0.026/bin/ifort -vec-report0
#FC = /opt/intel/bin/ifort
#FC = /opt/pgi/linux86-64/8.0-4/bin/pgf90 # c101
FC = /usr/local/bin/mpif90
#FC = gfortran

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
#INCLUDE = /usr/local/lib

INCLUDE = /usr/local/include
# Compiler flags
DBUGFLAGS = -g -traceback -check # debug version
FCFLAGS = -traceback -O3 -I$(INCLUDE) # release version
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

SOBJECTS = $(OTHERS) $(SERIALFN).o init_serial.o

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

release: $(POBJECTS)
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

preserial: 
	cp mendel.f90 $(SERIALFN).f90
	cat diagnostics.f90 >> $(SERIALFN).f90
	cp init.f90 init_serial.f90
	sed -i -e '/START_MPI/,/END_MPI/d' init_serial.f90
	sed -i -e '/START_MPI/,/END_MPI/d' $(SERIALFN).f90

cln:
	\rm -f mendel.o migration.o mendel

clean:
	\rm -f *.o *.mod $(TARGET) test mendel_serial.f90 migration_serial.f90 *.f90-e a.out\
	       success

###########################################
# dependencies
###########################################

sort.o:		sort.f90
	$(FC) $(FCFLAGS) -c sort.f90

random_pkg.o:	random_pkg.f90
	$(FC) $(FCFLAGS) -c random_pkg.f90

mendel.o:       mendel.f90 common.h
	$(FC) $(FCFLAGS) -c mendel.f90

init.o:		init.f90 common.h
	$(FC) $(FCFLAGS) -c init.f90

init_serial.o:		init_serial.f90 common.h
	$(FC) $(FCFLAGS) -c init_serial.f90

diagnostics.o:  diagnostics.f90 common.h
	$(FC) $(FCFLAGS) -c diagnostics.f90

fileio.o:	fileio.f90 common.h
	$(FC) $(FCFLAGS) -c fileio.f90

selection.o:    selection.f90 common.h
	$(FC) $(FCFLAGS) -c selection.f90

mating.o:	mating.f90
	$(FC) $(FCFLAGS) -c mating.f90

mutation.o:	mutation.f90
	$(FC) $(FCFLAGS) -c mutation.f90

migration.o:	migration.f90 common.h
	$(FC) $(FCFLAGS) -c migration.f90

mendel_serial.o:	mendel_serial.f90 common.h
	$(FC) $(FCFLAGS) -c mendel_serial.f90

migration_serial.o:	migration_serial.f90 common.h
	$(FC) $(FCFLAGS) -c migration_serial.f90

polygenic.o:	polygenic.f90
	$(FC) $(FCFLAGS) -c polygenic.f90

test.o:		test.f90
	$(FC) $(FCFLAGS) -c test.f90

profile.o:	profile.f90
	$(FC) $(FCFLAGS) -c profile.f90

inputs.o:       inputs.f90
	        $(FC) $(FCFLAGS) -c inputs.f90

genome.o:       genome.f90
	        $(FC) $(FCFLAGS) -c genome.f90


