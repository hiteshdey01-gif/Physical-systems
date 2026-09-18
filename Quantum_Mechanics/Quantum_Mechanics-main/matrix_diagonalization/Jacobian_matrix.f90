program jacobi_eigen
implicit none

integer, parameter :: n = 5
integer :: i,j,p,q,iter,max_iter

real(8) :: A(n,n), V(n,n)
real(8) :: d(n), b(n), z(n)
real(8) :: sm

real(8) :: t,c,s,tau
real(8) :: theta
real(8) :: g,h

max_iter = 100

! Read matrix from file


open(unit=10,file='hamiltonian_matrix.dat',status='old')

do i=1,n
    read(10,*) (A(i,j), j=1,n)
end do

close(10)


! Initialize eigenvector matrix


V = 0.0d0

do i=1,n
    V(i,i) = 1.0d0
end do

! Initialize diagonal arrays


do i=1,n
    b(i) = A(i,i)
    d(i) = A(i,i)
    z(i) = 0.0d0
end do

! Jacobi iteration


do iter = 1, max_iter

    sm = 0.0d0

    do p = 1, n-1
        do q = p+1, n
            sm = sm + abs(A(p,q))
        end do
    end do

    if (sm == 0.0d0) exit

    do p = 1, n-1
        do q = p+1, n

            if (abs(A(p,q)) > 0.0d0) then

                h = d(q) - d(p)

                if (abs(h) + abs(A(p,q)) == abs(h)) then
                    t = A(p,q)/h
                else
                    theta = 0.5d0*h/A(p,q)
                    t = 1.0d0/(abs(theta)+sqrt(1.0d0+theta*theta))
                    if (theta < 0.0d0) t = -t
                end if

                c = 1.0d0/sqrt(1.0d0+t*t)
                s = t*c
                tau = s/(1.0d0+c)

                h = t*A(p,q)

                z(p) = z(p) - h
                z(q) = z(q) + h

                d(p) = d(p) - h
                d(q) = d(q) + h

                A(p,q) = 0.0d0

                ! j < p
                do j = 1, p-1
                    g = A(j,p)
                    h = A(j,q)
                    A(j,p) = g - s*(h + g*tau)
                    A(j,q) = h + s*(g - h*tau)
                end do

                ! p < j < q
                do j = p+1, q-1
                    g = A(p,j)
                    h = A(j,q)
                    A(p,j) = g - s*(h + g*tau)
                    A(j,q) = h + s*(g - h*tau)
                end do

                ! j > q
                do j = q+1, n
                    g = A(p,j)
                    h = A(q,j)
                    A(p,j) = g - s*(h + g*tau)
                    A(q,j) = h + s*(g - h*tau)
                end do

                ! update eigenvectors
                do j = 1, n
                    g = V(j,p)
                    h = V(j,q)
                    V(j,p) = g - s*(h + g*tau)
                    V(j,q) = h + s*(g - h*tau)
                end do

            end if

        end do
    end do

    ! sweep correction
    do i = 1, n
        b(i) = b(i) + z(i)
        d(i) = b(i)
        z(i) = 0.0d0
    end do

end do


! Print Results


print *, "Eigenvalues:"
do i = 1, n
    print *, d(i)
end do

print *
print *, "Eigenvectors (columns correspond to eigenvalues):"

do i = 1, n
    write(*,'(5f12.6)') (V(i,j), j=1,n)
end do

end program jacobi_eigen  



