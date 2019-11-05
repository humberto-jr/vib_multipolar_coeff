real(kind = 8) function pes(sml_r, big_r, theta)
    implicit none

!   NOTE: on entry, sml_r and big_r have the same units used in the
!	wavefunctions, and theta has units of rad. Units of energy used
!	herein are also, by default, the units for all results.

    real(kind = 8), intent(in) :: sml_r, big_r, theta

    stop "pes(): you have to define the PES in pes.f90"

    pes = 1.0d0 ! Setup your PES here and remove the 'stop' above
end function pes
