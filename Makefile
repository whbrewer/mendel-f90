SRC = src
TARGET = $(SRC)/mendel
INSTALL_DIR = /usr/local/bin

FC = mpif90
INCLUDE ?= /usr/local/include
GIT_TAG := $(shell git describe --tags --abbrev=0 2>/dev/null)
GIT_SHA := $(shell git rev-parse --short=4 HEAD 2>/dev/null)
GIT_VERSION := $(if $(GIT_TAG),$(patsubst v%,%,$(GIT_TAG))-$(GIT_SHA),$(GIT_SHA))

# gfortran flags (LEGACYFLAGS needed for older Fortran compatibility)
FCFLAGS = -O3 -I$(SRC) -I$(INCLUDE) -J$(SRC)
DBUGFLAGS = -g -fbacktrace -fcheck=all -Wall

# Object files
MODULES = mpi_helpers sort random_pkg inputs genome profile polygenic init selection
CORE = $(MODULES) mutation mating fileio
OBJS = $(addprefix $(SRC)/, $(addsuffix .o, $(CORE) diagnostics mendel migration))
TEST_OBJS = $(addprefix $(SRC)/, $(addsuffix .o, $(CORE) diagnostics test migration))

.PHONY: all release debug test install uninstall dist clean cln

all: release

release: pre-build $(OBJS)
	$(FC) $(FCFLAGS) -o $(TARGET) $(OBJS)

debug: FCFLAGS = $(DBUGFLAGS)
debug: pre-build $(OBJS)
	$(FC) $(FCFLAGS) -o $(TARGET) $(OBJS)

test: pre-build $(TEST_OBJS)
	$(FC) $(FCFLAGS) -o test $(TEST_OBJS)

pre-build:
	@printf 'character(len=64), parameter :: build_version = "%s"\n' "$(GIT_VERSION)" > $(SRC)/version.inc
	@sed "s/tagsinput('add', '[^']*')/tagsinput('add', '$(GIT_VERSION)')/" spc-app/mendel.j2 > spc-app/mendel.j2.tmp
	@mv spc-app/mendel.j2.tmp spc-app/mendel.j2

install:
	install $(TARGET) $(INSTALL_DIR)

uninstall:
	rm -f $(INSTALL_DIR)/mendel

dist:
	scripts/package.sh

cln:
	rm -f $(SRC)/mendel.o $(SRC)/migration.o $(TARGET)

clean:
	rm -f $(SRC)/*.o $(SRC)/*.mod $(SRC)/version.inc $(TARGET) test a.out success
	$(MAKE) -C tests clean

# Pattern rule for Fortran compilation
$(SRC)/%.o: $(SRC)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

# Dependencies on common.h and mpi_helpers module
$(SRC)/mendel.o: $(SRC)/common.h $(SRC)/mpi_helpers.o
$(SRC)/init.o: $(SRC)/common.h $(SRC)/mpi_helpers.o
$(SRC)/diagnostics.o: $(SRC)/common.h $(SRC)/mpi_helpers.o
$(SRC)/migration.o: $(SRC)/common.h $(SRC)/mpi_helpers.o
$(SRC)/fileio.o: $(SRC)/common.h
$(SRC)/selection.o: $(SRC)/common.h
$(SRC)/mpi_helpers.o: $(SRC)/common.h
