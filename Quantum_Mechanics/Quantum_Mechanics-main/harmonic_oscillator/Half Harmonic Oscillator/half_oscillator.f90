! 1D Quantum half harmonic oscillator
program half_harmonic_oscillator_quantum
    
    implicit real*8(a-h,o-z)
    character ch
    common ak,al,h,n
    v(x)=0.5*ak*x**2  !potential energy
    pi=acos(-1.0)

    ak=0.5
    al=10.0
    h=0.0001
    eps=1.0e-4
    E=0.01

! Initial guess
    E3=3.0
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
    E1=E2
    E2=E3
    goto 5
    endif
10   continue
write(*,*)Anorm(E3)

call wavefunction(E3)
end program half_harmonic_oscillator_quantum




function Del(E)
    implicit real*8(a-h,o-z)
    common ak,al,h,n
   
    call right(E,siR,siRh)

    Del = si0
    return
end function Del

    
subroutine right(E,siR,siRh)
    implicit real*8(a-h,o-z)    
    common ak,al,h,n
    v(x)=0.5*ak*x**2
    x=al
    si0=0.0
    si1=1.0
    xt= sqrt(2*E/ak)
    nr= int((al/h +.5))
    do i=1,nr+1
        si2= (2.0- h*h*(E-v(x)))*si1 -si0
        si0=si1
        si1=si2
        x=x-h
    enddo
   
    si0=si1
    return
    end subroutine right

!Normalization of the wavefunction
function Anorm(E)
    implicit real*8(a-h,o-z)    
    common ak,al,h,n
    v(x)=0.5*ak*x**2
        
    nr= int(al/h+0.5)
    x=al
    si0=0.0 
    si1=1.0
    do i=1,nr+1
        si2= (2.0- h*h*(E-v(x)))*si1 -si0
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
    common ak,al,h,n
    v(x)=0.5*ak*x**2
   
   
   
    open(1,file='half.dat',status='unknown')
    A=Anorm(E)
   
    
    
    nr= int((al-sqrt(2*E/ak))/h+0.5)
    x=al
    si0=0.0
    si1=1.0
    do i=1,nr+1
        write(1,*) x, A*s*si0
        si2= (2.0- h*h*(E-v(x)))*si1 -si0
        si0=si1
        si1=si2
        x=x-h
    enddo
    close(1)
    return
    end subroutine wavefunction    


