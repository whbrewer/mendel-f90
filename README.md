# README #

Mendel's Accountant (MENDEL) is an advanced numerical simulation program for modeling genetic change over time and was developed collaboratively by Sanford, Baumgardner, Brewer, Gibson and ReMine.

For more information visit http://www.mendelsaccountant.info or http://sourceforge.net/projects/mendelsaccount

### Installation ###

* Build instructions: (1) edit Makefile, (2) run "make"

* Dependencies:  

  - Fortran compiler is required for compiling.  A free gfortran compile can easily by installed on a Linux system by executing the command "sudo apt-get install gfortran" on Debian/Ubuntu systems, and "sudo yum install gfortran" on Redhat/Centos/SuSe systems.  gfortran can also easily on Mac OS X using Homebrew (see brew.sh) by running "brew install gcc" (gfortran is included with gcc).

  - MPICH (optional). Since mendel uses mpich libraries, for parallel computations, you may download and install mpich from www.mpich.org/downloads/.  

  - However, it is possible to install mendel without installing the parallel libs, which means all options will work, but will be limited to running a single tribe/deme at a time. In order to compile Mendel without the parallel libs, first remove all MPICH code by running "make preserial", then run "make serial".  Sometimes the last compilation step needs to be run manually