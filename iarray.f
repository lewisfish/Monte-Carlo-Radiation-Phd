      subroutine iarray(xface,yface,zface,rhokap)

      include 'grid.txt'

      integer i,j,k

c**** Initialize array values to be zero

      do i=1,nxg+1
        xface(i)=0.
      end do
      do i=1,nyg+1
        yface(i)=0.
      end do
      do i=1,nzg+1
        zface(i)=0.
      end do

      do i=1,nxg
         do j=1,nyg
            do k=1,nzg
               rhokap(i,j,k)=0.
            end do
          end do
      end do

      return
      end
