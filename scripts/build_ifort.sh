#!/bin/bash

# modules
ifort -O3 -fpp -c legendre_polynomial_lib.f90
ifort -O3 -fpp -c quadpack_double_lib.f90
ifort -O3 -fpp -c vib_wavefunction.f90

# stand alone routines
ifort -O3 -fpp -c driver.f90 -fopenmp
ifort -O3 -fpp -c rot_integral.f90
ifort -O3 -fpp -c open_outputs.f90
ifort -O3 -fpp -c print_output.f90
ifort -O3 -fpp -c uvip3p.f90
ifort -O3 -fpp -c main.f90
ifort -O3 -fpp -c pes.f90

# executable
ifort -o main.out *.o -fopenmp
