INCLUDE ?= /usr/local/include
SRC_DIR := ../src
MODDIR := $(SRC_DIR)
PROJECT_INCLUDE := $(SRC_DIR)
DBUGFLAGS ?= -g -traceback -check $(LEGACYFLAGS)
FCFLAGS ?= -O3 -I$(PROJECT_INCLUDE) -I$(INCLUDE) -I$(MODDIR) -J$(MODDIR)
FC = mpif90

TARGET := test_main

SRC_OBJS := mpi_helpers.o inputs.o profile.o genome.o random_pkg.o polygenic.o \
	init.o unittest.o sort.o mutation.o selection.o migration.o fileio.o
TEST_OBJS := $(TARGET).o
OBJECTS := $(SRC_OBJS) $(TEST_OBJS)
TEST_PROGRAMS := test_inputs test_poisson test_polygenic test_restart test_competition

.PHONY: all clean tests run

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $(OBJECTS)

vpath %.f90 $(SRC_DIR)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

tests: $(TEST_PROGRAMS)

test_inputs: test_inputs.o inputs.o unittest.o
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $^

test_poisson: test_poisson.o random_pkg.o
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $^

test_polygenic: test_polygenic.o
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $^

test_restart: test_restart.o $(SRC_OBJS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $^

test_competition: test_competition.o $(SRC_OBJS)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $^

run: all tests
	./$(TARGET)
	./test_inputs
	./test_poisson
	./test_polygenic
	./test_restart
	./test_competition

clean:
	rm -f *.o *.mod $(TARGET) $(TEST_PROGRAMS) *.dmp.*
