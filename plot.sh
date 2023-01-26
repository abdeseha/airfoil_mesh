#!/bin/bash

gfortran mesh.f90
./a.out

gnuplot <<'EOF'
plot "Data.dat" w l
set term png
set output "airfoil.png"
replot
pause 20
EOF





