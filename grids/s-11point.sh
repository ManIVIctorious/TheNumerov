#!/bin/bash




for i in grid_*.dat;
    do
      ~/numerov_2d_multipoint_arpack-sparsefill/numerov2d_11p  $i 1.0  
    done



