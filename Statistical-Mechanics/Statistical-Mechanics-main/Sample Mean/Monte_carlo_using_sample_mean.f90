program monte carlo
implicit real*8(a-h,o-z)
f(x) = sin(x)**2

open(1,file='samplemean.dat')
pi= acos(-1.0)
a=0.0
b=2*pi

N = 1000
do j=1,6
  sum = 0.0
  do i= 1,N
    x= a+(b-a)*rand()
    sum = sum+f(x)
  enddo
result = (b-a)*sum/N
write(*,*)'the output is:',result
error = abs(pi-result)
write(1,*)log(float(N)),log(error)
N = N*10
enddo
stop
end program monte carlo
