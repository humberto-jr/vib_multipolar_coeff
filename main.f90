!
! Written by Humberto Jr (humberto.dasilvajr@unlv.edu)
! Oct, 2019
!

program main
    use vib_wavefunction, only: load_vib_wavef, free_resources
    implicit none

    logical :: use_omp
    real(kind = 8) :: R_min, R_max, R_step

    integer :: grid_size, n

    namelist /grid/ R_min, R_max, R_step, use_omp
    read(5, grid)

    if (R_min .gt. R_max) then
        stop "main(): R_min > R_max (input error)"
    end if

    grid_size = int((R_max - R_min)/R_step)

    if (grid_size < 1) then
        stop "main(): (R_max - R_min)/R_step < 1 (input error)"
    end if

    call load_vib_wavef()

    do n = 1, grid_size
        call driver(R_min + dble(n - 1)*R_step, use_omp)
    end do

    call free_resources()
end program
