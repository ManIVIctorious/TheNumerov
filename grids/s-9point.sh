#!/bin/bash




for i in grid_*.dat;
    do
      ~/numerov_2d_multipoint_arpack-sparsefill/numerov2d_09p  $i 1.0  
    done



