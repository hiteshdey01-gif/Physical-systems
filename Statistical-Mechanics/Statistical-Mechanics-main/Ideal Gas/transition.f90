program monte_carlo
    implicit none

    
    integer, parameter :: N=10, tmax=20
    integer :: ntrial, trial, t, i, index, NL 
    integer :: state(N)
    real*8 :: avg_NL(0:tmax)

    ntrial = 100000

    ! Initialize averages 
    do t = 0, tmax
        avg_NL(t) = 0.0
    end do  

    ! monte carlo simulation
    do trial = 1, ntrial
        do i = 1, N
            state(i) = 1     ! Initial state: all particles in the left box
        end do

        NL = N  ! Number of particles in the left box
        avg_NL(0) = avg_NL(0) + NL

        ! Time evolution
        do t = 1, tmax
            index = int(N * rand()) + 1  ! Randomly select a particle
            if (state(index).eq.1) then
                state(index) = 0  ! Move particle to the right box
            else
                state(index) = 1  ! Move particle to the left box
            end if

            NL = 0
            do i = 1,N
                NL = NL + state(i)  ! Count number of particles in the left box
            end do
            avg_NL(t) = avg_NL(t) + NL
        end do
    end do

    ! Average over trials
    do t = 0, tmax
        avg_NL(t) = avg_NL(t) / ntrial
    end do

    ! Output results
    open(1, file='avg_NL.dat')
    do t = 0, tmax
        write(1,*) t, avg_NL(t)
    end do
    close(1)
    write(*,*) 'Simulation completed. Results saved in avg_NL.dat'
    stop
end program monte_carlo
