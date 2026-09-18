      PROGRAM ELECTRONIC_LENS_AUTO

      IMPLICIT REAL*8(A-H,O-Z)


! --- GRID ---


      DIMENSION Phi(0:50,-40:40)
      DIMENSION X(0:5), Y(0:5)


! --- PARAMETERS ---


      m = 10
      n = 20
      mmax = 50
      nmax = 40

      Phi0 = 1000.0
      w    = 1.5
      eps  = 1.0E-7 * Phi0


! --- BOUNDARY CONDITIONS ---


      do k = 0, nmax
         Phi(mmax,k) = Phi0
      enddo

      do i = m, mmax-1
         Phi(i,n) = Phi0
      enddo

      do k = n, nmax
         Phi(m,k) = Phi0
      enddo

      do i = 0, mmax
         Phi(i,nmax) = i * (Phi0/m)
      enddo


! --- INITIAL GUESS ---


      do i = m, mmax-1
         do k = 0, n-1
            Phi(i,k) = Phi0
         enddo
      enddo

      do i = 0, m-1
         do k = 0, nmax-1
            Phi(i,k) = i * (Phi0/m)
         enddo
      enddo


! --- SOR ITERATION ---


      itn = 0

10   epsmax = 0.0

      do k = 0, nmax-1

         if (k .LT. n) then
            imax = mmax-1
         else
            imax = m-1
         endif

         do i = 1, imax

            if (k .EQ. 0) then
               U = (Phi(i+1,k) + Phi(i-1,k) + 4.0*Phi(i,k+1)) / 6.0
            else
               U = (Phi(i,k+1) + Phi(i,k-1) +  Phi(i+1,k) + Phi(i-1,k)) / 4.0+ (Phi(i,k+1) - Phi(i,k-1)) / (8.0*k)
            endif

            Phiold = Phi(i,k)
            Phi(i,k) = Phiold + w*(U - Phiold)

            err = ABS(Phi(i,k) - Phiold)
            if (err .GT. epsmax) epsmax = err

         enddo
      enddo

      itn = itn + 1

      if (epsmax .GT. eps) GOTO 10

      write(*,*) "Converged in iterations =", itn





      do k = -nmax, -1
         do i = 0, mmax
            Phi(i,k) = Phi(i,-k)
         enddo
      enddo


! AUTOMATIC EQUIPOTENTIAL GENERATION ---


      h = 1.0

    ! Loop over desired potentials (i = 2,4,6,8)
      do ichoice = 2, 8, 2

         p = Phi(ichoice,nmax)

    ! Create filename automatically
         if (ichoice .EQ. 2) OPEN(11, FILE='v200.dat', STATUS='UNKNOWN')
         if (ichoice .EQ. 4)  OPEN(11, FILE='v400.dat',STATUS='UNKNOWN')
         if (ichoice .EQ. 6) OPEN(11, FILE='v600.dat', STATUS='UNKNOWN')
         if (ichoice .EQ. 8) OPEN(11, FILE='v800.dat', STATUS='UNKNOWN')

         
       
        
        
    ! Top point
         z = ichoice*h
         r = nmax*h

         write(11,*)  z,  r
         write(11,*) -z,  r
         write(11,*)  z, -r
         write(11,*) -z, -r

! Loop over rows
         do k = nmax-1, 0, -1

            do i = 0, 4
               X(i) = Phi(ichoice+i,k)
               Y(i) = (ichoice+i)*h
            enddo

            z = Pn(X,Y,p)
            r = k*h

            write(11,*)  z,  r
            write(11,*) -z,  r
            write(11,*)  z, -r
            write(11,*) -z, -r

         enddo

         CLOSE(11)

         write(*,*) "Saved:", fname

      enddo

      STOP
      END


    ! --- INTERPOLATION FUNCTION ---


      FUNCTION Pn(X,Y,p)

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(0:5), Y(0:5), D(0:5,0:5)

      N = 4

      CALL Difftab(X,Y,D)

      sum = D(0,0)
      prod = 1.0

      do i = 1, N-1
         prod = prod*(p - X(i-1))
         sum = sum + prod*D(i,i)
      enddo

      Pn = sum
      RETURN
      END


! --- DifFERENCE TABLE ---


      SUBROUTINE Difftab(X,Y,D)

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(0:5), Y(0:5), D(0:5,0:5)

      N = 4

      do i = 0, N-1
         D(0,i) = Y(i)
      enddo

      do i = 1, N-1
         do j = i, N-1
            D(i,j) = (D(i-1,j) - D(i-1,j-1)) /(X(j) - X(j-i))
         enddo
      enddo

      RETURN
      END program 
