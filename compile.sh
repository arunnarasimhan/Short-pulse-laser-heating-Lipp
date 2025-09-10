#! /bin/bash

gfortran -ffree-line-length-0 -O3 -ffast-math -fexpensive-optimizations  -o nTTM_Si.x nTTM_Si.f90

## For debugging use the following instead:

#gfortran \
#	\
#	-ffree-line-length-0 -Wextra -Wall -pedantic -fbounds-check -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid,denormal \
#	-Wno-tabs -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused \
#	\
#	-o nTTM_Si.x nTTM_Si.f90
 
