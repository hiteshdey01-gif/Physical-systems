program hit_or_miss

    implicit real*8(a-h,o-z)

    N = 1000
    pi = acos(-1.0)

    open(1,file='hit_or_miss.dat')

    do j = 1,6
        
        hits = 0

        do i = 1,N
            call random_number(r1)
            call random_number(r2)
            x = -2.0 + 4.0*r1
            y = -1.0 + 2.0*r2
            if ((x*x)/4.0 + y*y .le.1.0) then
                hits = hits + 1
            endif
        enddo

        area = 8.0*hits/N
        exact_area = 2.0 *pi
        error = abs(area - exact_area)

        write(*,*) 'Estimated area = ', area
        write(1,*) log(float(N)),log(error)
        N = N*10
    enddo
    close(1)
    stop
end program hit_or_miss
