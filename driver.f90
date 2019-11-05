subroutine driver(big_r, use_omp)
    use vib_wavefunction, only: v_min, v_max
    use omp_lib
    implicit none

    logical, intent(in) :: use_omp
    real(kind = 8), intent(in) :: big_r

    integer :: n, v, v_prime
    real(kind = 8) :: start_time, end_time

    logical, save :: first_call = .true., is_homo = .false.
    integer, save :: lambda_max = 0, lambda_step = 1, n_max

    integer, allocatable, save :: index_a(:), index_b(:)

    namelist /rotor/ lambda_max, is_homo

    if (first_call) then
        first_call = .false.

        rewind(5)
        read(5, rotor)

        if (is_homo) then
            call open_outputs(v_min, v_max, lambda_max, 2)
        else
            call open_outputs(v_min, v_max, lambda_max, 1)
        end if

!       NOTE: number of upper diagonal elements: n*(n + 1)/2
        n_max = (v_max - v_min + 1)*(v_max - v_min + 2)/2

        allocate(index_a(1:n_max), index_b(1:n_max))

        if (.not. allocated(index_a) .or. .not. allocated(index_b)) then
            stop "driver(): unable to allocate resources"
        end if

        n = 0
        do v = v_min, v_max
            do v_prime = v, v_max
                n = n + 1
                index_a(n) = v
                index_b(n) = v_prime
            end do
        end do

        if (n .ne. n_max) then
            stop "driver(): n .ne. n_max (internal error)"
        end if
    end if

    !$omp parallel do default(none) shared(index_a, index_b, lambda_max, is_homo, big_r, n_max) private(start_time, end_time) schedule(static) if(use_omp)
        do n = 1, n_max
            start_time = omp_get_wtime()
            call print_output(index_a(n), index_b(n), lambda_max, is_homo, big_r)
            end_time = omp_get_wtime()

            !$omp critical
                write(6, fmt = *) "#"
                write(6, fmt = *) "# v      = ", index_a(n)
                write(6, fmt = *) "# v'     = ", index_b(n)
                write(6, fmt = *) "# R      = ", big_r
                write(6, fmt = *) "# wtime  = ", end_time - start_time, " s"
                write(6, fmt = *) "# thread = ", omp_get_thread_num()
            !$omp end critical
        end do
    !$omp end parallel do
end subroutine driver
