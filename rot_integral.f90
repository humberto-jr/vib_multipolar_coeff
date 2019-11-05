subroutine rot_integral(a, b, lambda_max, is_homo, big_r, result)
    use quadpack_double_lib, only: dqag
    implicit none

    logical, intent(in) :: is_homo
    real(kind = 8), intent(in) :: big_r
    integer, intent(in) :: a, b, lambda_max

    real(kind = 8), intent(out) :: result(0:lambda_max)

    real(kind = 8) :: theta_max, error, scale_fact
    integer :: max_step, info, last, lambda, lambda_step

    integer, parameter :: limit = 10000
    integer, parameter :: lenw = 4*limit
    real(kind = 8), parameter :: pi = acos(-1.0d0)

    integer, allocatable :: iwork(:)
    real(kind = 8), allocatable :: work(:), legendre_poly(:, :)

    allocate(iwork(limit))

    if (.not. allocated(iwork)) then
        stop "rot_integral(): unable to allocate resources (1)"
    end if

    allocate(work(lenw))

    if (.not. allocated(work)) then
        stop "rot_integral(): unable to allocate resources (2)"
    end if

    allocate(legendre_poly(1, 0:lambda_max))

    if (.not. allocated(legendre_poly)) then
        stop "rot_integral(): unable to allocate resources (3)"
    end if

    if (is_homo) then
        lambda_step = 2
        scale_fact  = 2.0d0
        theta_max   = pi/2.0d0
    else
        lambda_step = 1
        scale_fact  = 1.0d0
        theta_max   = pi
    end if

    result = 0.0d0

    do lambda = 0, lambda_max, lambda_step

!       NOTE: gauss-kronrod integration rule with 30-61 points
        call dqag(f, 0.0d0, theta_max, 0.1d0, 1.0d-12, 5, result(lambda), error, max_step, info, limit, lenw, last, iwork, work)

        result(lambda) = scale_fact*dble(2*lambda + 1)*result(lambda)/2.0d0

        if (info .ne. 0) then
            write(6, fmt = *) "rot_integral(): dqag() failed with error code = ", info
            stop
        end if
    end do

    deallocate(iwork, work, legendre_poly)
    return

    contains
    real(kind = 8) function f(x)
        use vib_wavefunction, only: qag_vib_quadrature
        use legendre_polynomial_lib, only: p_polynomial_value
        implicit none

        real(kind = 8), intent(in) :: x

        real(kind = 8) :: vib_integral, cosx(1)

!       NOTE: dqag() will invoke f(x) with x (theta) in rad
        call qag_vib_quadrature(a, b, big_r, x, vib_integral)

        cosx(1) = cos(x)
        call p_polynomial_value(1, lambda_max, cosx, legendre_poly)

        f = vib_integral*legendre_poly(1, lambda)*sin(x)
    end function f

end subroutine rot_integral
