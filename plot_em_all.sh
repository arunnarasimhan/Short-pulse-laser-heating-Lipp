#! /bin/bash

{ gnuplot plot_Ne_contour.gpi;
convert -density 300 +antialias contour_Ne.eps contour_Ne.png;
rm contour_Ne.eps;} &


{ gnuplot plot_Te_contour.gpi; 
convert -density 300 +antialias contour_Te.eps contour_Te.png; 
rm contour_Te.eps; } &


{ gnuplot plot_Ta_contour.gpi;
convert -density 300 +antialias contour_Ta.eps contour_Ta.png;
rm contour_Ta.eps; } &


{ gnuplot plot_energy.gpi;
convert -density 300 +antialias Ea-and-Ee.eps Ea-and-Ee.png;
rm Ea-and-Ee.eps; } &


{ gnuplot plot_error.gpi;
convert -density 300 +antialias Err.eps Err.png;
rm Err.eps; } 
