#!/bin/bash

# modules
gfortran -O3 -ffree-line-length-none -c legendre_polynomial_lib.f90
gfortran -O3 -ffree-line-length-none -c quadpack_double_lib.f90
gfortran -O3 -ffree-line-length-none -c vib_wavefunction.f90

# stand alone routines
gfortran -O3 -ffree-line-length-none -c driver.f90 -fopenmp
gfortran -O3 -ffree-line-length-none -c rot_integral.f90
gfortran -O3 -ffree-line-length-none -c open_outputs.f90
gfortran -O3 -ffree-line-length-none -c print_output.f90
gfortran -O3 -ffree-line-length-none -c uvip3p.f90
gfortran -O3 -ffree-line-length-none -c main.f90
gfortran -O3 -ffree-line-length-none -c pes.f90

# executable
gfortran -o main.out *.o -fopenmp
