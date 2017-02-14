#!/bin/bash




for i in morse_*.dat;
    do
      ~/numerov1d_multipoint_gls/numerov_gsl_14 $i 9.411764705882e-01 13 
    done



