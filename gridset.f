      subroutine gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,
     +                  kappa1,kappa2,id)

      implicit none

      include 'grid.txt'

      real xmax,ymax,zmax,kappa1,kappa2

      integer i,j,k,kflag,id
      real x,y,z,rho,taueq1,taupole1,taueq2,taupole2
      if(id.eq.0.)then
      print*, ' '
      print *, 'Setting up density grid....'
      end if
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
           call density(x,y,z,rho,kflag,kappa1,kappa2,i,j,k,rhokap)
           rhokap(i,j,k,1)=rho
           kflag=0
           call density(x,y,z,rho,kflag,kappa1,kappa2,i,j,k,rhokap)
           rhokap(i,j,k,2)=rho                    
        end do
       end do
      end do

c****************** Calculate equatorial and polar optical depths ****
      taueq1=0.
      taupole1=0.
      taueq2=0.
      taupole2=0.
      do i=1,nxg
         taueq1=taueq1+rhokap(i,nyg/2,nzg/2,1)
         taueq2=taueq2+rhokap(i,nyg/2,nzg/2,2)
      enddo
      do i=1,nzg
         taupole1=taupole1+rhokap(nxg/2,nyg/2,i,1)
         taupole2=taupole2+rhokap(nxg/2,nyg/2,i,2)
      enddo
      taueq1=taueq1*2.*xmax/nxg
      taupole1=taupole1*2.*zmax/nzg
      taueq2=taueq2*2.*xmax/nxg
      taupole2=taupole2*2.*zmax/nzg
      if(id.eq.0.)then
      print *,'taueq1 = ',taueq1,'  taupole1 = ',taupole1
      print *,'taueq2 = ',taueq2,'  taupole2 = ',taupole2
      end if
      open(10,file='density1.dat')
      open(11,file='density2.dat')
      open(12,file='density3.dat')
      do i=1,nxg
           write(10,*) (rhokap(i,102,j,1),j=1,nzg)
           write(11,*) (rhokap(i,j,102,1),j=1,nzg)
           write(12,*) (rhokap(102,i,j,1),j=1,nzg)
           end do
      close(10)
      close(11)
      

      return
      end
