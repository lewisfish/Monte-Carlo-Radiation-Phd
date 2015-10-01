      subroutine gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,
     +                  kappa1,kappa2)

      implicit none

      include 'grid.txt'

      real xmax,ymax,zmax,kappa1,kappa2

      integer i,j,k,kflag
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
           kflag=1
           call density(x,y,z,rho,kflag,kappa1,kappa2)
           rhokap(i,j,k,1)=rho
           kflag=0
           call density(x,y,z,rho,kflag,kappa1,kappa2)
           rhokap(i,j,k,2)=rho                    
        end do
       end do
      end do

c****************** Calculate equatorial and polar optical depths ****
!      taueq=0.
!      taupole=0.
!      do i=1,nxg
!         taueq=taueq+rhokap(i,nyg/2,nzg/2)
!      enddo
!      do i=1,nzg
!         taupole=taupole+rhokap(nxg/2,nyg/2,i)
!      enddo
!      taueq=taueq*2.*xmax/nxg
!      taupole=taupole*2.*zmax/nzg
!      print *,'taueq = ',taueq,'  taupole = ',taupole

      open(10,file='density1.dat')
      open(11,file='density2.dat')
      do i=1,nxg
           write(10,*) (rhokap(i,100,j,1),j=1,nzg)
           write(11,*) (rhokap(i,100,j,2),j=1,nzg)
           end do
      close(10)
      close(11)
      

      return
      end
