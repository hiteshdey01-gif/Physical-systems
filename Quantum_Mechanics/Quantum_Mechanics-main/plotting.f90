program plotting_ground_wavefunction
    implicit none

    integer,parameter :: Nd=5
    integer :: i

    real(8) :: E0
    real(8) :: xmin,xmax,x,h
    real(8) :: psi_total   ! final wavefunction
    real(8) :: psi  
    integer :: k, npoints
    

    real(8), dimension(Nd) :: C  ! Eigenvector coefficients

    open(10, file='jacobi_eigen.dat', status='old')

    read(10,*) E0
    do i=1,Nd
        read(10,*) C(i)
    enddo
    close(10)
    print *, "Ground state energy: ", E0
    print *, "Eigenvector coefficients: ", C(1:Nd)

    ! Define the range and step size for plotting
    xmin = -8.0
    xmax = 8.0
    h = 0.05

    open(20, file='ground_wavefunction.dat', status='unknown')
    
    npoints = (xmax - xmin)/h
    do k=0,npoints
        x = xmin + k*h
        psi_total = 0.0
        do i = 1,Nd
            psi_total = psi_total + C(i)*psi(i-1,x)  ! i-1 because wavefunction index starts from 0
        enddo
        write(20,*) x, psi_total
    enddo
    close(20)
end program plotting_ground_wavefunction


real(8) function factorial(n)
    implicit none
    integer :: n,i
    factorial = 1.0
    do i=1,n
        factorial = factorial * i
    enddo
    end function factorial

    real(8) function hermite(n,x)
    implicit none
    integer :: n,i
    real(8) :: x
    real(8) :: h0,h1,h2
    if (n.eq.0) then
        hermite = 1.0
        return
    endif
    if (n.eq.1) then
        hermite = 2.0 * x
        return
    endif
    h0 = 1.0
    h1 = 2.0 * x    
    do i=2,n
        h2 = 2.0*x*h1 - 2.0*(i-1)*h0
        h0 = h1
        h1 = h2 
    enddo
    hermite = h1
    end function hermite

    real(8) function psi(n,x)
    implicit none
    integer :: n
    real(8) :: x
    real(8) :: norm
    real(8) :: hermite, factorial
    norm = 1.0 / sqrt((2.0**n)*factorial(n)*sqrt(acos(-1.0)))
    psi = norm * hermite(n,x) * exp(-0.5*x**2)
    end function psi





