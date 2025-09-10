# Short-pulse-laser-heating-Lipp

Instructions and general guidelines for running this code.

1. Open the "nTTM_Si.f90" file and execute "./compile.sh"
2. One helpful way to store the data for easy collection and plotting is to type the following "./nTTM_SI.x > case3ChenAndBeraun.out". This will place the code output into a file named "case3ChenAndBeraun.out". The program stores the history of carrier density (n), electron temperature (Te) and lattice temperature Tl in 3 columns. Having this data in text file format makes it easy to retrieve and plot.


Make sure to have a fortran compiler installed (gfortran). Other useful VS code extensions for fortran are "Code Runner" and "Modern Fortran". These make it very easy to run fortran code. 
