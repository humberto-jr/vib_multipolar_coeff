subroutine open_outputs(v_min, v_max, lambda_max, lambda_step)
    implicit none

    integer, intent(in) :: v_min, v_max, lambda_max, lambda_step

    character(len = 10) :: a, b, c
    character(len = 1024) :: filename

    integer :: v, v_prime, lambda, file_exists

    do v = v_min, v_max
        write(a, fmt = '(i5)') v
        a = adjustl(a)

        do v_prime = v, v_max
            write(b, fmt = '(i5)') v_prime
            b = adjustl(b)

            do lambda = 0, lambda_max, lambda_step
                write(c, fmt = '(i5)') lambda
                c = adjustl(c)

                filename = "v="//trim(a)//"-"//trim(b)//"_lambda="//trim(c)//".dat"

!               NOTE: delete old files, if any
                open(unit = 666, file = trim(filename), iostat = file_exists, status = "old")
                if (file_exists .eq. 0) close(666, status = "delete")

!               NOTE: open brand new output files
                open(unit = 666, file = trim(filename), action = "write", status = "new")
                close(666)
            end do
        end do
    end do
end subroutine open_outputs
