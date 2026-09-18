program harmonic_oscillator
    implicit none
    integer :: n, i, No  ! n- quantum number, i- loop counter, No- number of points
    real(8) :: x, xmin, xmax, dx, psi
    real(8) :: wavefunction  ! declare functions

    open(10, file='hermite.dat', status='unknown')  ! Open file to write results

    n=4 ! Choose quantum number
    xmin = -5.0
    xmax = 5.0
    No = 500  ! for smoth curve we choose 500 points
    dx = (xmax - xmin) / (No )

    do i= 0,No  ! we start from 0 so that x=xmin
        x= xmin + i*dx
        psi = wavefunction(n,x)  ! calculate wavefunction value
        write(10,*) x, psi  ! write x and psi to file
    enddo
    close(10)  ! Close the file
end program harmonic_oscillator


    ! Subprogram to calculate Hermite polynomial
    real(8) function hermite(n,x)
    implicit none
    integer :: n,k ! order of polynomial, loop index
    real(8) :: x   ! coordinate
    real(8) :: h0,h1,h2   ! H(n-1), H(n), H(n+1)- temporary value

    ! case 1 n=0
    if (n.eq.0) then
        hermite = 1.0
        return
    endif

    ! case 2 n=1
    if (n.eq.1) then
        hermite = 2.0 * x
        return  
    endif

    ! if n>=2 we use recursion relation. we first asssign the two polynomials
    h0 = 1.0       ! H(0)
    h1 = 2.0 * x   ! H(1)
    do k=2,n
        h2 = 2.0*x*h1 - 2.0*(k-1)*h0  ! recursion relation
        h0 = h1  ! update H(n-1)
        h1 = h2  ! update H(n)
    enddo
    hermite = h1  ! return H(n)
    end function hermite

    real (8) function factorial(n)
        implicit none
        integer :: n, i ! n- factorial order, i- loop counter
        real(8) :: f
        f = 1.0
        do i=1,n
            f = f * i
        enddo
        factorial = f   ! return factorial value
    end function factorial
    ! even if n=0 the code  will work correctly, 
    !because the do loop will not execute and f will remain 1.0, which is the correct value for 0!

    ! Wavefunction Function
    real (8) function wavefunction(n,x)
        implicit none
        integer :: n
        real(8) :: x
        real(8) :: Hn, fact, norm
        real(8) :: hermite, factorial  ! declare functions
        real(8), parameter :: pi = 3.14159265358979323846
        Hn = hermite(n,x)  ! calculate Hermite polynomial
        fact = factorial(n)  ! calculate n!
        norm = (1.0/sqrt((2.0**n)*fact))*(1.0/pi)**0.25
        
        wavefunction = norm * Hn *exp(-0.5*x**2)  ! return wavefunction value
    end function wavefunction
