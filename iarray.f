      subroutine iarray(xface,yface,zface,rhokap)

      include 'grid.txt'

      integer i,j,k,p

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
              do p=1,2
               rhokap(i,j,k,p)=0.
               end do
            end do
          end do
      end do

      return
      end
