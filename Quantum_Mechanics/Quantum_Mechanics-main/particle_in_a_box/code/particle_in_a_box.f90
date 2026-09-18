! one dimensionnal particle in a box
program pbox
    implicit real*8(a-h,o-z)
    character ch
    common h,n

    pi=acos(-1.0)
    aL=10.0 ! box length
    h= 0.01 ! step size
    epsilon = 1.0e-4 ! convergencee criteria
    n=int(aL/h)+1 ! number of points

!    Secant method 
    E3= 0.0
3    E1 = E3 + 0.1 ! initial guess
    E2 = E1 + 0.1 ! next initial guess

5    f1 = sil(E1)
    f2 = sil(E2)
    E3 = (f1*E2 - f2*E1)/(f1 - f2) ! next guess
    check = abs(E3 - E2)

    if (check.le.epsilon) then
        write(*,*) 'Eigenvalue converged: ', E3
        goto 10
    else
        E1 = E2
        E2 = E3
        goto 5
    end if

10   continue
     call wavefunction(E3)
    ! write(*,*)'Plot wavefunction.dat'
   ! write(*,*)'Want to get next excited state? (y/n)'
    !read(*,*) ch
    !if (ch .eq. 'y' .or. ch .eq. 'Y') then
    !    E3 = E3 + 0.1 ! next initial guess for the next excited state
    !    goto 3 
    !endif


     stop
end program pbox


! Shooting method
function sil(E)
    implicit real*8(a-h,o-z)
    common h,n

    x = 0.0
    si0 = 0.0 ! initial condition
    si1 = 1.0 ! arbitrarily chosen
    do i= 1,n
        si2 =(2-h**2*E)*si1 - si0 ! finite difference formula
        si0 = si1
        si1 = si2
        x = x + h
    end do
    sil = si1
    return
end function sil

! Normalization of the wavefunction
function Anorm(E)
    implicit real*8(a-h,o-z)
    common h,n

    x = 0.0
    si0 = 0.0 ! initial  boundary condition
    si1 = 1.0 ! arbitrarily chosen
    sum = 0.0

    do i= 1,n
        si2 =(2-h**2*E)*si1 - si0 ! finite difference formula
        si0 = si1
        si1 = si2
        x = x + h
        sum = sum + si1**2
    end do

    Anorm = 1.0/sqrt(sum*h) ! normalization constant
    return
end function Anorm

! Normalized wavefunction
subroutine wavefunction(E)
    implicit real*8(a-h,o-z)
    common h,n

    x = 0.0
    si0 = 0.0 ! initial  boundary condition
    si1 = 1.0 ! arbitrarily chosen
    open(1, file='wavefunction.dat', status='unknown') ! open file to save wavefunction
    An = Anorm(E) ! normalization constant

    do i= 1,n
        write(1,*) x, An*si0 ! save normalized wavefunction to file
        si2 =(2-h**2*E)*si1 - si0 ! finite difference formula
        si0 = si1
        si1 = si2
        x = x + h
       
    end do
    close(1)
    return

end subroutine wavefunction

       
