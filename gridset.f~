      subroutine gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,kappa)

      implicit none

      include 'grid.txt'

      real xmax,ymax,zmax,kappa

      integer i,j,k
      real x,y,z,rho,taueq,taupole

      print *, 'Setting up density grid....'

c**********  Linear Cartesian grid. Set up grid faces ****************
      do i=1,nxg+1
         xface(i)=(i-1)*2.*xmax/nxg
      end do
      do i=1,nyg+1
         yface(i)=(i-1)*2.*ymax/nyg
      end do
      do i=1,nzg+1
         zface(i)=(i-1)*2.*zmax/nzg
      end do

c**************  Loop through x, y, and z to set up grid density.  ****
      do i=1,nxg
       do j=1,nyg
        do k=1,nzg
           x=xface(i)-xmax+xmax/nxg
           y=yface(j)-ymax+ymax/nyg
           z=zface(k)-zmax+zmax/nzg

c**********************Call density setup subroutine 
           call density(x,y,z,rho)
           rhokap(i,j,k)=rho*kappa

        end do
       end do
      end do

c****************** Calculate equatorial and polar optical depths ****
      taueq=0.
      taupole=0.
      do i=1,nxg
         taueq=taueq+rhokap(i,nyg/2,nzg/2)
      enddo
      do i=1,nzg
         taupole=taupole+rhokap(nxg/2,nyg/2,i)
      enddo
      taueq=taueq*2.*xmax/nxg
      taupole=taupole*2.*zmax/nzg
      print *,'taueq = ',taueq,'  taupole = ',taupole
      
c************** Write out density grid as unformatted array
      open(10,file='density.dat',status='unknown',
     +        form='unformatted')
           write(10) rhokap
      close(10)

      return
      end
