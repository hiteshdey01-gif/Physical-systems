program gaussian_wave_packet

    implicit real*8(a-h,o-z)

    integer :: nx,nk,nt
    integer :: i,j,m
    real*8 :: xmin,xmax,dx
    real*8 :: kmin,kmax,dk
    real*8 :: tmax,dt
    real*8 :: x,k,t
    real*8 :: a, k0
    real*8 :: phi, omega
    real*8 :: psir, psii
    real*8 :: sumr, sumi
    real*8 :: weight
    real*8 :: prob
    

    parameter (nx=401)
    parameter (nk=801)
    parameter (nt=80)

    dimension xr(nx)
    

    a= 4.0
    k0 = 1.5
  
    xmin = -40.0
    xmax = 40.0
    dx = (xmax-xmin)/(nx-1)
    do i= 1,nx
        xr(i) = xmin+(i-1)*dx
    enddo
    
    
    kmin = -10.0
    kmax = 10.0
    dk = (kmax-kmin)/(nk-1)


    tmax = 30.0
    dt = tmax/(nt-1)

    open(10,file='gaussian_wave_packet.dat')


    do m=1,nt
        t=(m-1)*dt
        do i=1,nx
            x=xr(i)
            sumr = 0.0
            sumi = 0.0

            do j= 1,nk
                k = kmin+(j-1)*dk
                phi = exp(-(a*a*(k-k0)**2)/4.0) 
                omega = 0.5*k**2

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

            sumr = (sumr*dk/3.0)
            sumi = (sumi*dk/3.0)

            prob = sumr*sumr + sumi*sumi
            write(10,*) x,t,prob
        enddo
        write(10,*)
        
    enddo
    close(10)

end program gaussian_wave_packet

