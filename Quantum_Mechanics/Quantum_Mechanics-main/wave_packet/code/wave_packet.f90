program wavepacket
    implicit none

    ! Parameters
    ! Number of spatial points, momentum points, and time steps
    integer, parameter :: Nx = 400
    integer, parameter :: Nk = 600
    integer, parameter :: Nt = 50

    real(8), parameter :: a = 1.0   !  wave packet width

    ! spatial range and momentum range
    real(8), parameter :: xmin = -10.0
    real(8), parameter :: xmax = 10.0
    real(8), parameter :: kmin = -20.0
    real(8), parameter :: kmax = 20.0

    real(8) :: dx,dk

    ! Defining arrays to store spatial and momentum grid points
    real(8) :: x(Nx), k(Nk)

    real(8) :: pi
    complex(16) :: psi(Nx)
    real(8) :: phi(Nk)
    integer :: i, j
    real(8) :: prob   ! probability density

     real(8) :: tmax, t, dt  ! time parameters
    real(8) :: omega         ! angular frequency for time evolution
    complex(16) :: sum,term   ! temporary variables for summation in time evolution
    integer :: n             

    ! Initialize the spatial and momentum grids
    dx = (xmax - xmin) / (Nx - 1) !  using (Nx-1) ensures that the last point equals xmax. So the grid is perfectly symmetric.
    dk = (kmax - kmin) / (Nk - 1)

    ! Initialize the spatial grid
    do i = 1, Nx
        x(i) = xmin + (i - 1) * dx ! if we h used x(i) = xmin + i*dx, then x(1) would be xmin + dx, which is not what we want. We want x(1) to be exactly xmin.
    end do

    ! Initialize the momentum grid
    do j = 1, Nk
        k(j) = kmin + (j - 1) * dk
    end do

   
    pi = acos(-1.0)

    ! compute phi(k) 
    ! special case for k=0 to avoid division by zero.
    do j = 1,Nk
        if (abs(k(j)) .lt. 0.00001) then
            phi(j) = sqrt(2.0/(pi*a)) * (a/2)
        else
            phi(j) = sqrt(2.0/(pi*a)) * sin(k(j)*a/2)/k(j)
        endif
    end do
    ! So now the program has computed the initial momentum space wave function phi(k) for a rectangular wave packet of width a.
    ! The special case for k=0 ensures that we get the correct limit of phi(k) as k approaches zero.
    ! Physically this corresponds to the fact that a rectangular wave packet has a sinc-shaped momentum distribution.


    ! using nested loops
        ! loop over time
            ! loop over position
                ! compute psi(x,t) using the k- integraton


    ! time parameters
    tmax = 5.0
    dt = tmax / Nt

    open(10, file='wavepacket.dat', status='unknown') ! open a file to write the results


    do n= 1,Nt ! loop over time steps
        t = (n-1) * dt
        do i = 1, Nx  ! loop over spatial positions
            sum = (0.0, 0.0)  ! sum is a complex number, first part is real, second part is imaginary.
            do j = 1,Nk  ! integration over momentum k
                omega = k(j)**2 / 2.0
                term = phi(j) * exp((0.0, 1.0) * (k(j)*x(i) - omega*t))
                if ( j ==1 .or. j == Nk) then  ! Simpson integration method: first and last terms
                    sum = sum + term
                else if (mod(j,2).eq.0) then  ! mod(j,2) ==0 checks if the number is even
                    sum = sum + 4.0*term
                else    ! odd terms
                    sum = sum + 2.0*term
                end if
                
            end do
            psi(i) = (dk/3.0)*sum/sqrt(2.0*pi)

            ! compute probability density
            prob = real(psi(i))**2 + aimag(psi(i))**2  ! |psi|^2 = Re(psi)^2 + Im(psi)^2
            write(10,*) x(i), prob
        end do
        write(10,*)     ! blank line to separate time steps in the output file
        write(10,*)     ! In gnuplot we can use "index" to separate different time steps when plotting.And two consecutive blank lines seprate data blocks.each block corresopnd to one time step.
    end do
    close(10)
end program wavepacket

    
