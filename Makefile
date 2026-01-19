INSTALL_DIR = /usr/local/bin
SRC_DIR = src
OBJDIR = $(SRC_DIR)
MODDIR = $(SRC_DIR)
PROJECT_INCLUDE = $(SRC_DIR)

GIT_VERSION := $(shell git describe --abbrev=7 --dirty --always --tags)

FC = mpif90

# Optional MPI libs (usually provided by mpif90).
#LIBS = -lmpich -lpthread -lmpl
#LIBS = -L/usr/local/lib -lmpich -lpthread -lmpl
# Open MPI
#LIBS = -L/usr/lib64/mpi/gcc/openmpi/lib64 -lopen-rte -lmpi

INCLUDE ?= /usr/local/include
# Compiler flags
LEGACYFLAGS ?= -fallow-argument-mismatch -std=legacy
DBUGFLAGS = -g -traceback -check $(LEGACYFLAGS) # debug version
#FCFLAGS = -traceback -O3 -I$(INCLUDE) # release version ifort
FCFLAGS = -O3 -I$(PROJECT_INCLUDE) -I$(INCLUDE) -I$(MODDIR) -J$(MODDIR) $(LEGACYFLAGS) # release version gfortran
# note use flag -fpe:0 to handle floating point exceptions

# Linker flags
#LDFLAGS = -static-libgfortran -static-libgcc

# executable name
TARGET = $(SRC_DIR)/mendel

MODULES = $(OBJDIR)/mpi_helpers.o $(OBJDIR)/sort.o $(OBJDIR)/random_pkg.o \
          $(OBJDIR)/inputs.o $(OBJDIR)/genome.o $(OBJDIR)/profile.o \
          $(OBJDIR)/polygenic.o $(OBJDIR)/init.o $(OBJDIR)/selection.o

OTHERS = $(MODULES) $(OBJDIR)/mutation.o $(OBJDIR)/mating.o \
         $(OBJDIR)/fileio.o

POBJECTS = $(OTHERS) $(OBJDIR)/diagnostics.o $(OBJDIR)/mendel.o \
           $(OBJDIR)/migration.o

TOBJECTS = $(OTHERS) $(OBJDIR)/diagnostics.o $(OBJDIR)/test.o \
           $(OBJDIR)/migration.o

##########################################
# build rules
##########################################

.PHONY: all debug test release parallel install uninstall dist clean cln

all: release

debug: FCFLAGS = $(DBUGFLAGS)
debug: pre-build $(POBJECTS)
	$(FC) $(DBUGFLAGS) $(LDFLAGS) -o $(TARGET) $(POBJECTS) $(LIBS)

test: pre-build $(TOBJECTS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o test $(TOBJECTS) $(LIBS)

pre-build:
	printf 'character(len=64), parameter :: build_version = "%s"\n' "$(GIT_VERSION)" > $(SRC_DIR)/version.inc

release: pre-build $(POBJECTS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $(TARGET) $(POBJECTS) $(LIBS)

parallel: release

install:
	install $(TARGET) $(INSTALL_DIR)

uninstall:
	echo "removing file $(INSTALL_DIR)/$(TARGET)"
	rm $(INSTALL_DIR)/$(TARGET)

dist:
	scripts/package.sh

cln:
	\rm -f $(OBJDIR)/mendel.o $(OBJDIR)/migration.o $(TARGET)

clean:
	\rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(SRC_DIR)/version.inc $(TARGET) test0* *.f90-e a.out success
	$(MAKE) -C tests clean

###########################################
# dependencies
###########################################

$(OBJDIR)/sort.o:		$(SRC_DIR)/sort.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/sort.f90 -o $(OBJDIR)/sort.o

$(OBJDIR)/random_pkg.o:	$(SRC_DIR)/random_pkg.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/random_pkg.f90 -o $(OBJDIR)/random_pkg.o

$(OBJDIR)/mendel.o:       $(SRC_DIR)/mendel.f90 $(PROJECT_INCLUDE)/common.h $(OBJDIR)/mpi_helpers.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/mendel.f90 -o $(OBJDIR)/mendel.o

$(OBJDIR)/init.o:		$(SRC_DIR)/init.f90 $(PROJECT_INCLUDE)/common.h $(OBJDIR)/mpi_helpers.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/init.f90 -o $(OBJDIR)/init.o

$(OBJDIR)/diagnostics.o:  $(SRC_DIR)/diagnostics.f90 $(PROJECT_INCLUDE)/common.h $(OBJDIR)/mpi_helpers.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/diagnostics.f90 -o $(OBJDIR)/diagnostics.o

$(OBJDIR)/fileio.o:	$(SRC_DIR)/fileio.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/fileio.f90 -o $(OBJDIR)/fileio.o

$(OBJDIR)/selection.o:    $(SRC_DIR)/selection.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/selection.f90 -o $(OBJDIR)/selection.o

$(OBJDIR)/mating.o:	$(SRC_DIR)/mating.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/mating.f90 -o $(OBJDIR)/mating.o

$(OBJDIR)/mutation.o:	$(SRC_DIR)/mutation.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/mutation.f90 -o $(OBJDIR)/mutation.o

$(OBJDIR)/migration.o:	$(SRC_DIR)/migration.f90 $(PROJECT_INCLUDE)/common.h $(OBJDIR)/mpi_helpers.o
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/migration.f90 -o $(OBJDIR)/migration.o

$(OBJDIR)/polygenic.o:	$(SRC_DIR)/polygenic.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/polygenic.f90 -o $(OBJDIR)/polygenic.o

$(OBJDIR)/test.o:		$(SRC_DIR)/test.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/test.f90 -o $(OBJDIR)/test.o

$(OBJDIR)/profile.o:	$(SRC_DIR)/profile.f90
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/profile.f90 -o $(OBJDIR)/profile.o

$(OBJDIR)/inputs.o:       $(SRC_DIR)/inputs.f90
	        $(FC) $(FCFLAGS) -c $(SRC_DIR)/inputs.f90 -o $(OBJDIR)/inputs.o

$(OBJDIR)/genome.o:       $(SRC_DIR)/genome.f90
	        $(FC) $(FCFLAGS) -c $(SRC_DIR)/genome.f90 -o $(OBJDIR)/genome.o
$(OBJDIR)/mpi_helpers.o:	$(SRC_DIR)/mpi_helpers.f90 $(PROJECT_INCLUDE)/common.h
	$(FC) $(FCFLAGS) -c $(SRC_DIR)/mpi_helpers.f90 -o $(OBJDIR)/mpi_helpers.o
