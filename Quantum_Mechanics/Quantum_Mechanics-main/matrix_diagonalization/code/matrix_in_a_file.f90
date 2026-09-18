program matrix
    dimension A(10,10)
    open(1,file='matrix.dat')

    do i=1,5
        do j=1,5
            A(i,j) = i*j
        enddo
    write(1,*)(int(A(i,j)),j=1,5)
    enddo
    close(1)
end program matrix        
