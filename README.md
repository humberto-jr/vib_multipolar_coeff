AUTHOR
======

Written by Humberto Jr (humberto.dasilvajr@unlv.edu), Oct, 2019.

Based on the original routine vlambda.f90 from Fabio Carelli (2015-2016).



About
=====

This program computes the vibrational multipolar coefficients as function of the
scattering coordinate, R, for an atom-diatom problem, defined as

V(v', v, lambda, R) := [(2*lambda + 1)/2]*
integral[V(v', v, R, theta)*P(lambda, cos(theta))*sin(theta)*dtheta] from 0 to pi

where,

V(v', v, R, theta) := <v'|V(r, R, theta)|v>

{|v>} = set of vibrational wavefunctions of the diatomic target

V(r, R, theta) = potential energy surface of the system

P(lambda, x) = Lengendre polynomial of order lambda evaluated at x



HOW TO USE
==========

1) Vibrational wavefunctions are expected in independent data files named
v=0.dat, v=1.dat, v=2.dat etc, each of which with two columns of data: r and the
probability amplitude. All blank lines or lines starting with '#' are ignored by
the code.

2) The potential energy surface (PES) in Jacobi coordinates (r, R, theta) is
defined by the user in pes.f90. On entry, units for (r, R) are the same of those
used in the wavefunctions. Theta is in rad. On exit, the energy unit chosen is
also the unit for all results.

3) Invoke either build_gfortran.sh or build_ifort.sh to compile the program,
e.g.

./build_gfortran.sh

4) Write an input file, e.g.

&grid
    R_min   = 0.5
    R_max   = 10.0
    R_step  = 0.25
    use_omp = .true.
&end

&vib_wavefunction
    v_min = 0
    v_max = 1
	v_dir = "/home/my_name/my_calculations/my_wavefunctions"
&end

&rotor
    lambda_max = 4
    is_homo = .false.
&end

Where, 'use_omp' switch on/off the use of OpenMP; 'v_dir' is the folder in which
the wavefuctions are stored; and, 'is_homo' tells the program to consider a
(hetero) homonuclear diatomic target.

5) Invoke the program, e.g.

./main.out < my_input_file.d

6) Each multipolar coefficient, as function of R, will be printed in independent
files named v=0-0_lambda=0.dat, v=0-1_lambda=0.dat, etc. In which units are
driven by the wavefunctions and the pes (see above).

7) One is required to write his own post-processing script to read the output
data files and rewrite in the format used by ASPIN or any other scattering code.
