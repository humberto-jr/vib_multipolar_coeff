subroutine print_output(a, b, lambda_max, is_homo, big_r)

    logical, intent(in) :: is_homo
    real(kind = 8), intent(in) :: big_r
    integer, intent(in) :: a, b, lambda_max

    integer :: lambda, lambda_step

    character(len = 1024) :: filename
    character(len = 10) :: a_char, b_char, c_char

    real(kind = 8), allocatable :: result(:)

    allocate(result(0:lambda_max))

    if (.not. allocated(result)) then
        stop "print_output(): unable to allocate resources"
    else
        result = 0.0d0
    end if

    call rot_integral(a, b, lambda_max, is_homo, big_r, result)

    if (is_homo) then
        lambda_step = 2
    else
        lambda_step = 1
    end if

    write(a_char, fmt = '(i5)') a
    a_char = adjustl(a_char)

    write(b_char, fmt = '(i5)') b
    b_char = adjustl(b_char)

    do lambda = 0, lambda_max, lambda_step
        write(c_char, fmt = '(i5)') lambda
        c_char = adjustl(c_char)

        filename = "v="//trim(a_char)//"-"//trim(b_char)//"_lambda="//trim(c_char)//".dat"

        open(unit = 666, file = trim(filename), access = "append", action = "write", status = "old")
        write(666, fmt = '(F20.8, 4x, ES20.12)') big_r, result(lambda)
        close(666)
    end do

    deallocate(result)
end subroutine print_output
