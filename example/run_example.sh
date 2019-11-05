#!/bin/bash

# These are the steps to run the example of He + H2+. The input file (He+H2+_example.d),
# potential energy surface (pes.f90) and the H2+ wavefunctions (H2+_vib_wavefunctions/)
# are pre-built.

# Step 1: to compile the code (only object files, no executable at this point)
cd ../
./build_gfortran.sh
mv *.o *.out *.mod example/
cd example/

# Step 2: to compile our own pes.f90 for He + H2+
gfortran -O3 -ffree-line-length-none -c pes.f90

# Step 3: to link the new pes.o object file to the rest of the code
gfortran -o main.out *.o -fopenmp
rm *.o *.mod

# Step 4: run the program
./main.out < He+H2+_example.d > He+H2+_example.log

echo "All done. See He+H2+_example.log and all generated .dat files."
