Author
------

Written by Humberto Jr (humberto.dasilvajr@unlv.edu).

Based on the original routine vlambda.f90 from Fabio Carelli (2015-2016).

Oct, 2019


About
-----

This program computes the vibrational multipolar coefficients as function of the
scattering coordinate, R, for an atom-diatom problem, defined as

![alt text](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Logo Title Text 1")

where,

1. The phi functions are rovibrational wavefunctions with quantum number upsilon,
from the diatomic target, as function of the internuclear distance r.

2. The V function is a potential energy surface of the system in Jacobi coordinates (r, R, theta).

3. And P represents Lengendre polynomials of order lambda evaluated at cos(theta).



How to use
----------

1. Vibrational wavefunctions are expected in independent data files named
v=0.dat, v=1.dat, v=2.dat etc, each of which with two columns of data: r and the
probability amplitude. All blank lines or lines starting with '#' are ignored by
the code.

2. The potential energy surface (PES) in Jacobi coordinates (r, R, theta) is
defined by the user in pes.f90. On entry, units for (r, R) are the same of those
used in the wavefunctions. Theta is in rad. On exit, the energy unit chosen is
also the unit for all results.

3. Invoke either build_gfortran.sh or build_ifort.sh to compile the program, e.g.
./build_gfortran.sh (for GNU gfortran).

4. Write an input file, e.g. example/He+H2+_example.d; where, 'use_omp' switchs
on/off the use of OpenMP; 'v_dir' is the folder in which the wavefuctions are
stored; and, 'is_homo' tells the program to consider a (hetero) homonuclear
diatomic target.

5. Invoke the program, e.g. ./main.out < my_input_file.d

6. Each multipolar coefficient, as function of R, will be printed in independent
files named v=0-0_lambda=0.dat, v=0-1_lambda=0.dat, etc. In which units are
driven by the wavefunctions and the PES (see above).

7. One is required to write his own post-processing script to read the output
data files and rewrite in the format used by ASPIN, Molscat or any other
scattering code.
