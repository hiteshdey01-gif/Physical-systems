program gravity_wave_packet

    implicit real*8(a-h,o-z)

    integer :: nx,nk,nt
    integer :: i,j,m
    real*8 :: xmin,xmax,dx
    real*8 :: kmin,kmax,dk
    real*8 :: tmax,dt
    real*8 :: x,k,t
    real*8 :: a,lambda0, k0, g
    real*8 :: phi, omega
    real*8 :: psir, psii
    real*8 :: sumr, sumi
    real*8 :: weight
    real*8 :: prob
    real*8 :: pi,norm

    parameter (nx=401)
    parameter (nk=801)
    parameter (nt=50)

    dimension xr(nx)
    pi = 4.0*atan(1.0)


    a= 6.0
    lambda0 = 10.0
    k0 = 2.0*pi/lambda0
    g=9.8
    norm = 1.0/sqrt(2.0*pi)

    xmin = -50.0
    xmax = 50.0
    dx = (xmax-xmin)/(nx-1)
    do i= 1,nx
        xr(i) = xmin+(i-1)*dx
    enddo
    
    
    kmin = -5.0
    kmax = 5.0
    dk = (kmax-kmin)/(nk-1)


    tmax = 20.0
    dt = tmax/(nt-1)

    open(10,file='gravity.dat')


    do m=1,nt
        t=(m-1)*dt
        do i=1,nx
            x=xr(i)
            sumr = 0.0
            sumi = 0.0

            do j= 1,nk
                k = kmin+(j-1)*dk
                phi = exp(-(a*a*(k-k0)**2)/4.0) + exp(-(a*a*(k+k0)**2)/4.0)
                omega = sqrt(g*abs(k))

                psir = phi*cos(k*x-omega*t)
                psii = phi*sin(k*x-omega*t)
                
                if(j.eq.1.or.j.eq.nk) then
                    weight= 1.0
                else if(mod(j,2).eq.0) then
                    weight = 4.0
                else
                    weight = 2.0
                endif
                
                sumr = sumr+ weight*psir
                sumi = sumi + weight*psii
            enddo

            sumr = (sumr*dk/3.0)*norm
            sumi = (sumi*dk/3.0)*norm

            prob = sumr*sumr + sumi*sumi
            write(10,*) x,t,prob
        enddo
        write(10,*)
        
    enddo
    close(10)

end program gravity_wave_packet

