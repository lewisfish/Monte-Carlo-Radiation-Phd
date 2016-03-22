MODULE gridset_mod

implicit none
save

CONTAINS
   subroutine gridset(id)

   use density_mod
   use constants, only : nxg,nyg,nzg,xmax,ymax,zmax
   use iarray, only    : rhokap,xface,yface,zface

   implicit none


   integer i,j,k,kflag,id
   DOUBLE PRECISION x,y,z,rho,taueq1,taupole1,taueq2,taupole2

   if(id.eq.0.)then
   print*, ' '
   print *, 'Setting up density grid....'
   end if
   !**********  Linear Cartesian grid. Set up grid faces ****************
   do i=1,nxg+1
      xface(i)=(i-1)*2.*xmax/nxg
   end do
   do i=1,nyg+1
      yface(i)=(i-1)*2.*ymax/nyg
   end do
   do i=1,nzg+1
      zface(i)=(i-1)*2.*zmax/nzg
   end do
    
   !**************  Loop through x, y, and z to set up grid density.  ****
   do i=1,nxg
    do j=1,nyg
     do k=1,nzg
        x=xface(i)-xmax+xmax/nxg
        y=yface(j)-ymax+ymax/nyg
        z=zface(k)-zmax+zmax/nzg
   !**********************Call density setup subroutine 
        kflag=1
        call density(x,y,z,rho)
        rhokap(i,j,k,kflag)=rho
     end do
    end do
   end do

   !****************** Calculate equatorial and polar optical depths ****
   taueq1=0.
   taupole1=0.
   taueq2=0.
   taupole2=0.
   do i=1,nxg
      taueq1=taueq1+rhokap(i,nyg/2,nzg/2,1)
   end do
   do i=1,nzg
      taupole1=taupole1+rhokap(nxg/2,nyg/2,i,1)
   end do
   taueq1=taueq1*2.*xmax/nxg
   taupole1=taupole1*2.*zmax/nzg
   if(id.eq.0.)then
   print *,'taueq1 = ',taueq1,'  taupole1 = ',taupole1
   end if

   end subroutine gridset
end MODULE gridset_mod
