program entropy_full

    implicit none

    integer :: n, No
    real(8) :: omega, P, S, total_states, binomial

    No = 10
    total_states = 2.0d0**No

    open(unit=10, file='entropy_full.dat', status='unknown')

    write(10,*) '# n        Omega(n)       P(n)       S(n)'

    do n = 0, No

        omega = binomial(No, n)
        P = omega / total_states
        S = log(omega)

        write(10,'(I3,3F15.6)') n, omega, P, S

    end do

    close(10)

    print*, 'Data written to entropy_full.dat'
end program entropy_full


!-------------------------------
! Factorial Function
!-------------------------------
    function factorial(k) result(fact)

        implicit none
        integer :: k, i
        real(8) :: fact

        fact = 1.0d0

        do i = 1, k
            fact = fact * dble(i)
        end do

    end function factorial

!-------------------------------
! Binomial Function
!-------------------------------
    function binomial(No, n) result(res)

        implicit none
        integer :: No, n
        real(8) :: res, factorial

        res = factorial(No) / ( factorial(n) * factorial(No - n) )

    end function binomial

