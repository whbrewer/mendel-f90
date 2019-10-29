# README #

**Note: this is the original/old version of Mendel's Accountant implemented in Fortran. This is not used anymore. It is kept here for reference.**

Mendel's Accountant (MENDEL) is an advanced numerical simulation program for modeling genetic change over time and was developed collaboratively by Sanford, Baumgardner, Brewer, Gibson and ReMine.

For more information visit http://www.mendelsaccountant.info or http://sourceforge.net/projects/mendelsaccount

### INSTALLATION ###

* Build instructions: (1) edit Makefile, (2) run "make"

* A Fortran compiler is required for compiling.  A free gfortran compile can easily by installed on a Linux system by executing the command "sudo apt-get install gfortran" on Debian/Ubuntu systems, and "sudo yum install gfortran" on Redhat/Centos/SuSe systems.  

* gfortran can also easily on Mac OS X using Homebrew (see brew.sh) by running "brew install gcc" (gfortran is included with gcc).  You'll also need to install two files: libmpfr.4.dylib and libmpc.3.dylib in /usr/local/lib.

* MPICH (optional). Since mendel uses mpich libraries, for parallel computations, you may download and install mpich from www.mpich.org/downloads/.  

* However, it is possible to install mendel without installing the parallel libs, which means all options will work, but will be limited to running a single tribe/deme at a time. In order to compile Mendel without the parallel libs, first remove all MPICH code by running "make preserial", then run "make serial".  It will croak when it tries to do the final linking.  So run this command after it fails:

gfortran -O3 -I/usr/local/include -static-libgfortran -static-libgcc -o mendel_serial sort.o random_pkg.o inputs.o genome.o profile.o polygenic.o selection.o mutation.o mating.o fileio.o  mendel_serial.o init_serial.o

### INTERFACE ###

Mendel's Accountant was designed to work with the Scientific Platform for the Cloud (SPC).  Instructions for setting up SPC can be found here: https://bitbucket.org/whbrewer/spc
