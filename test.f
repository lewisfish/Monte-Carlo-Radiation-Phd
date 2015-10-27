      program test

      implicit none

      integer :: i,j
      real,allocatable :: A

      allocate(A(1:2,1:2))

      do i=1,2
            do j=1,2
                  A(i,j)=A(i,j)+i
            end do
      end do

      print*, A

      end program test
