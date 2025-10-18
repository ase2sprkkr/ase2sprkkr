subroutine AAA(cdata, data, na, nb)
    complex*16 cdata(na, nb)
    real*8 data(na, nb)
    integer na
    integer nb
    integer i,j
    cdata = 2.0
    do i=1,na
        do j=1,nb
          data(i,j) = 1000*i+j
    end do
    end do

end subroutine
