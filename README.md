# Mendel's Accountant - Original Fortran Source Code

Mendel's Accountant (MENDEL) is an advanced numerical simulation program
for modeling genetic change over time. It was developed collaboratively
by Sanford, Baumgardner, Brewer, Gibson, and ReMine.

## Installation

- Build: edit `Makefile`, then run `make`.
- Linux (gfortran):
  - Debian/Ubuntu: `sudo apt-get install gfortran`
  - Redhat/CentOS/SUSE: `sudo yum install gfortran`
- macOS (Homebrew):
  - `brew install gcc` (includes `gfortran`)
  - `brew install open-mpi`
  - Install `libmpfr.4.dylib` and `libmpc.3.dylib` in `/usr/local/lib`
- MPICH (optional, for parallel runs): https://www.mpich.org/downloads/

## SPC integration

Mendel's Accountant was designed to work with the Scientific Platform for
the Cloud (SPC). Setup instructions:
https://github.com/whbrewer/spc

## Mendel's Accountant References

- Sanford et al., "Mendel's Accountant: A biologically realistic forward-time population genetics program." Scalable Computing: Practice and Experience 8, no. 2 (2007).
- Sanford et al., "Using computer simulation to understand mutation accumulation dynamics and genetic load." In International Conference on Computational Science, 386-392. Berlin, Heidelberg: Springer Berlin Heidelberg, 2007.
- Baumgardner et al., "Mendel’s Accountant: A new population genetics simulation tool for studying mutation and natural selection." In Proceedings of the Sixth International Conference on Creationism, vol. 8798. 2008.
- Sanford et al., "Using numerical simulation to test the validity of neo-Darwinian theory." In Proceedings of the International Conference on Creationism, vol. 6, no. 1, 16. 2008.
- Sanford et al., "Selection threshold severely constrains capture of beneficial mutations." In Biological Information: New Perspectives, 264-297. 2013.
- Brewer et al., "Information loss: potential for accelerating natural genetic attenuation of RNA viruses." In Biological Information: New Perspectives, 369-384. 2013.
- Brewer et al., "Using numerical simulation to test the 'mutation-count' hypothesis." In Biological Information: New Perspectives, 298-311. 2013.
- Gibson et al., "Can Purifying Natural Selection Preserve Biological Information?" In Biological Information: New Perspectives, 232-263. 2013.
- Baumgardner et al., "Can synergistic epistasis halt mutation accumulation? Results from numerical simulation." In Biological Information: New Perspectives, 312-337. 2013.
- Sanford et al., "The waiting time problem in a model hominin population." Theoretical Biology and Medical Modelling 12, no. 1 (2015): 18.
- Sanford et al., "Adam and Eve, designed diversity, and allele frequencies." 2018.
