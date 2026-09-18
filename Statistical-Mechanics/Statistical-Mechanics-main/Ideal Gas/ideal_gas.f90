program ideal_gas

    implicit real*8 (a-h,o-z)

    integer :: state(10)     ! array to hold the state of the system
    integer :: trial, ntrial, t
    
    ntrial = 100

    open(1, file='ideal_gas.dat', status='unknown')

    do j = 1, 6
        avg_NL = 0.0
        do trial = 1, ntrial
            ! Initialize all particles to LEFT
            do i = 1, 10
            state(i) = 1   ! 1 represents LEFT
            end do

          ! Repeat the process for a number of time unit to simulate the system

            do t = 1,10
                index = int(rand() * 10) + 1   
                if (state(index) .eq. 1) then
                    state(index) = 0   ! Flip to RIGHT
                else
                    state(index) = 1   ! Flip to LEFT
                end if
            end do

            ! Count the number of particles in LEFT and RIGHT
    
            count_left = 0
            do i = 1, 10
                count_left = count_left + state(i)   ! state(i) is 1 for LEFT and 0 for RIGHT
            end do

            avg_NL = avg_NL + count_left
        end do
        avg_NL = avg_NL / ntrial
        write(1,*) ntrial, avg_NL

        ntrial = ntrial * 10
    end do
    close(1)
    stop 

end program ideal_gas
