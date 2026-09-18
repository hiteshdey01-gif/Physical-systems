  ! 1D Quartic oscillator
program quartic_oscillator
    
    implicit real*8(a-h,o-z)
    character ch
    common ak,al,h,b
    v(x)=0.5*ak*x**2 +b*x**4
     
    b=0.2
    ak=0.5
    al=10.0
    h=0.001
    eps=1.0e-4
    E=0.01

! Initial guess
    E3=6.5
3    E1=E3+0.1
    E2=E1+0.1

! Secant Method
5   f1=Del(E1)
    f2=Del(E2)
    E3= (f1*E2-f2*E1)/(f1-f2)
    check = abs(E3-E2)
    if(check.le.eps)then
        write(*,*) 'Energy Eigenvalue =',E3
    goto 10
    else
    e1=E2
    E2=E3
    goto 5
    endif
10   continue
write(*,*)Anorm(E3)

call wavefunction(E3)
end program quartic_oscillator




function Del(E)
    implicit real*8(a-h,o-z)
    common ak,al,h,b

    call left(E,siL,siLh)
    call right(E,siR,siRh)

    s=siL/siR
    Del = siLh-s*siRh
    return
end function Del

subroutine left(E,siL,siLh)
    implicit real*8(a-h,o-z)
    common ak,al,h,b
    v(x)=0.5*ak*x**2 +b*x**4
    x=-al
    si0=0.0
    si1=1.0
    xt= sqrt(2*E/ak)
    nl= int((al+xt)/h+0.5)
    do i=1,nl+1
        ak0=2.0*(E-v(x-h))
        ak1=2.0*(E-v(x))
        ak2=2.0*(E-v(x+h))

        si2=(2.0*(1.0- 5*h*h*ak1/12.0)*si1 - (1.0 + h*h*ak0/12.0)*si0)/(1.0+h*h*ak2/12.0)
        si0=si1
        si1=si2
        x=x+h
    enddo
    siL=si0
    siLh=si1
    return
    end subroutine left
    
subroutine right(E,siR,siRh)
    implicit real*8(a-h,o-z)    
    common ak,al,h,b
    v(x)=0.5*ak*x**2 +b*x**4
    x=al
    si0=0.0
    si1=1.0
    xt= sqrt(2*E/ak)
    nr= int((al-xt)/h+0.5)
    do i=1,nr+1
        ak0=2.0*(E-v(x-h))
        ak1=2.0*(E-v(x))
        ak2=2.0*(E-v(x+h))

        si2=(2.0*(1.0- 5*h*h*ak1/12.0)*si1 - (1.0 + h*h*ak0/12.0)*si0)/(1.0+h*h*ak2/12.0)
        si0=si1
        si1=si2
        x=x-h
    enddo
    siR=si1
    siRh=si0
    return
    end subroutine right

!Normalization of the wavefunction
function Anorm(E)
    implicit real*8(a-h,o-z)    
    common ak,al,h,b
    v(x)=0.5*ak*x**2 +b*x**4
    x=-al
    si0=0.0
    si1=1.0
    sum=0.0
    xt= sqrt(2*E/ak)
    nl= int((al+sqrt(2*E/ak))/h+0.5)

    do i =1,nl+1
        ak0=2.0*(E-v(x-h))
        ak1=2.0*(E-v(x))
        ak2=2.0*(E-v(x+h))

        si2=(2.0*(1.0- 5*h*h*ak1/12.0)*si1 - (1.0 + h*h*ak0/12.0)*si0)/(1.0+h*h*ak2/12.0)
        si0=si1
        si1=si2
        x=x+h
        sum = sum + si1*si1
    enddo
    
    nr= int((al-sqrt(2*E/ak))/h+0.5)
    x=al
    si0=0.0 
    si1=1.0
    do i=1,nr+1
        ak0=2.0*(E-v(x-h))
        ak1=2.0*(E-v(x))
        ak2=2.0*(E-v(x+h))

        si2=(2.0*(1.0- 5*h*h*ak1/12.0)*si1 - (1.0 + h*h*ak0/12.0)*si0)/(1.0+h*h*ak2/12.0)
        si0=si1
        si1=si2
        x=x-h
        sum = sum + si1*si1
    enddo
    Anorm=1.0/sqrt(sum*h)
    return
    end function Anorm


subroutine wavefunction(E)
    implicit real*8(a-h,o-z)
    common ak,al,h,b
    v(x)=0.5*ak*x**2 +b*x**4
    call left(E,siL,siLh)
    call right(E,siR,siRh)
    s=siL/siR
    write(*,*)'The scale factor is',s
    open(1,file='oscq.dat',status='unknown')
    A=Anorm(E)
    
    
    x=-al
    si0=0.0
    si1=1.0
    nl= int((al+sqrt(2*E/ak))/h+0.5)
    xt= sqrt(2*E/ak)
    do i=1,nl+1
        write(1,*) x, A*si0
        ak0=2.0*(E-v(x-h))
        ak1=2.0*(E-v(x))
        ak2=2.0*(E-v(x+h))

        si2=(2.0*(1.0- 5*h*h*ak1/12.0)*si1 - (1.0 + h*h*ak0/12.0)*si0)/(1.0+h*h*ak2/12.0)
        si0=si1 
        si1=si2
        x=x+h   
    enddo
    
    nr= int((al-sqrt(2*E/ak))/h+0.5)
    x=al
    si0=0.0
    si1=1.0
    do i=1,nr+1
        write(1,*) x, A*s*si0
        ak0=2.0*(E-v(x-h))
        ak1=2.0*(E-v(x))
        ak2=2.0*(E-v(x+h))

        si2=(2.0*(1.0- 5*h*h*ak1/12.0)*si1 - (1.0 + h*h*ak0/12.0)*si0)/(1.0+h*h*ak2/12.0)
        si0=si1
        si1=si2
        x=x-h
    enddo
    close(1)
    return
    end subroutine wavefunction    

