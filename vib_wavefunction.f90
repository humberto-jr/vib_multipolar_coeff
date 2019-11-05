module vib_wavefunction
    implicit none

    public :: load_vib_wavef, qag_vib_quadrature, free_resources
    private :: init_vib_wavef, open_vib_wavef, vib_integrand

!   spline order
    integer, parameter, private :: order = 3

!   directory to search the wavefunctions
    character(len = 1024), private :: v_dir

!   range for all vib. quantum numbers, v (output)
    integer, public :: v_min, v_max

    type spline
        integer :: length
        real(kind = 8) :: x_min, x_max
        real(kind = 8), allocatable :: x(:), y(:)
    end type spline

!   spline interpolant
    class(spline), allocatable, private :: vib_wavef(:)

    contains
    subroutine init_vib_wavef()
        namelist /vib_wavefunction/ v_min, v_max, v_dir

        rewind(5)
        read(5, vib_wavefunction)

!       write(6, fmt = *) "# v (min, max) = ", v_min, v_max

        if (v_max < v_min) then
            stop "init_vib_wavef(): v_max < v_min (input error)"
        end if

        allocate(vib_wavef(v_min:v_max))

        if (.not. allocated(vib_wavef)) then
            stop "init_vib_wavef(): unable to allocate resources"
        end if
    end subroutine init_vib_wavef

!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!

    subroutine open_vib_wavef(input, v, max_row)

        integer, intent(in) :: input, v

        integer, intent(out) :: max_row

        character(len = 5) :: v_char
        character(len = 1024) :: filename, line

        write(v_char, '(i5)') v
        v_char = adjustl(v_char)
        filename = trim(v_dir)//"/v="//trim(v_char)//".dat"

        open(unit = input, file = trim(filename), action = "read", status = "old")

        max_row = 0
        do while (.true.)
            read(input, fmt = '(a)', end = 1) line

            if (line(1:1) .ne. '#' .and. len_trim(line) > 0) then
                max_row = max_row + 1
            end if
        end do

1       rewind(input)
        write(6, fmt = *) "# reading "//trim(filename)//"; num. of lines = ", max_row

        if (max_row < 1) then
            stop "open_vib_wavef(): invalid number of lines (input error)"
        end if
    end subroutine open_vib_wavef

!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!

    subroutine load_vib_wavef()

        integer :: v, n, n_max
        character(len = 1024) :: line

        call init_vib_wavef()

        do v = v_min, v_max
            call open_vib_wavef(666, v, n_max)

            allocate(vib_wavef(v)%x(1:n_max), vib_wavef(v)%y(1:n_max))

            if (.not. allocated(vib_wavef(v)%x) .or. .not. allocated(vib_wavef(v)%y)) then
                stop "load_vib_wavef(): unable to allocate resources"
            end if

            n = 0
            do while (n .ne. n_max)
                read(666, fmt = '(a)', end = 1) line

                if (line(1:1) .ne. '#' .and. len_trim(line) > 0) then
                    n = n + 1
                    read(line, fmt = *) vib_wavef(v)%x(n), vib_wavef(v)%y(n)
                end if
            end do

1           close(666)
            vib_wavef(v)%length = n_max
            vib_wavef(v)%x_min  = vib_wavef(v)%x(1)
            vib_wavef(v)%x_max  = vib_wavef(v)%x(n_max)
        end do
    end subroutine load_vib_wavef

!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!

    real(kind = 8) function vib_integrand(a, b, sml_r, big_r, theta)

        integer, intent(in) :: a, b
        real(kind = 8), intent(in) :: sml_r, big_r, theta

        real(kind = 8) :: x(1), wavef_a(1), wavef_b(1), pes, error

        x(1) = sml_r

        error = 0.0d0
        call uvip3p(order, vib_wavef(a)%length, vib_wavef(a)%x, vib_wavef(a)%y, 1, x, wavef_a, error)

        if (error .ne. 0.0d0) then
            write(6, fmt = *) "vib_integrand(): uvip3p() failed (1) with error code = ", error
            stop
        end if

        error = 0.0d0
        call uvip3p(order, vib_wavef(b)%length, vib_wavef(b)%x, vib_wavef(b)%y, 1, x, wavef_b, error)

        if (error .ne. 0.0d0) then
            write(6, fmt = *) "vib_integrand(): uvip3p() failed (2) with error code = ", error
            stop
        end if

        vib_integrand = wavef_a(1)*pes(sml_r, big_r, theta)*wavef_b(1)
    end function vib_integrand

!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!

    subroutine qag_vib_quadrature(a, b, big_r, theta, result)
        use quadpack_double_lib, only: dqag

        integer, intent(in) :: a, b
        real(kind = 8), intent(in) :: big_r, theta

        real(kind = 8), intent(out) :: result

        integer :: max_step, info, last
        real(kind = 8) :: r_min, r_max, error

        integer, parameter :: limit = 10000
        integer, parameter :: lenw = 4*limit

        integer, allocatable :: iwork(:)
        real(kind = 8), allocatable :: work(:)

        r_min = max(vib_wavef(a)%x_min, vib_wavef(b)%x_min)
        r_max = min(vib_wavef(a)%x_max, vib_wavef(b)%x_max)

        allocate(iwork(limit))

        if (.not. allocated(iwork)) then
            stop "qag_vib_quadrature(): unable to allocate resources (1)"
        end if

        allocate(work(lenw))

        if (.not. allocated(work)) then
            stop "qag_vib_quadrature(): unable to allocate resources (2)"
        end if

!       NOTE: gauss-kronrod integration rule with 30-61 points
        call dqag(f, r_min, r_max, 0.1d0, 1.0d-12, 5, result, error, max_step, info, limit, lenw, last, iwork, work)

        if (info .ne. 0) then
            write(6, fmt = *) "qag_vib_quadrature(): dqag() failed with error code = ", info
            stop
        end if

        deallocate(iwork, work)
        return

        contains
        real(kind = 8) function f(x)

            real(kind = 8), intent(in) :: x

            f = vib_integrand(a, b, x, big_r, theta)
        end function f

    end subroutine qag_vib_quadrature

!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!

    subroutine free_resources()

        deallocate(vib_wavef)
    end subroutine free_resources
end module vib_wavefunction
