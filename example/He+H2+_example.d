&grid
  R_min   = 0.9
  R_max   = 3.0
  R_step  = 0.025
  use_omp = .true.
&end

&vib_wavefunction
  v_min = 0
  v_max = 5
  v_dir = "./H2+_vib_wavefunctions"
&end

&rotor
  lambda_max = 6
  is_homo = .true.
&end
