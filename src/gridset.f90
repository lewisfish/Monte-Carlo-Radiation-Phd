MODULE gridset_mod

implicit none
save

CONTAINS
   subroutine gridset(id, wave)
   
   use ch_opt
   use density_mod
   use constants, only : nxg,nyg,nzg,xmax,ymax,zmax
   use iarray, only    : rhokap,xface,yface,zface,albedo,conc

   implicit none


   integer i,j,k,kflag,id
   DOUBLE PRECISION x,y,z,rho,taueq1,taupole1,taueq2,taupole2, wave

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
!   do i=1,nxg
!    do j=1,nyg
!     do k=1,nzg
!        x=xface(i)-xmax+xmax/nxg
!        y=yface(j)-ymax+ymax/nyg
!        z=zface(k)-zmax+zmax/nzg
!   !**********************Call density setup subroutine 
!        kflag=1
!        call density(x,y,z,rho)
!        rhokap(i,j,k,1)=rho
!     end do
!    end do
!   end do
   rhokap=0.d0
   albedo=0.d0
   do i=1,nzg
      z=zface(i)-zmax+zmax/nzg

      if(z.gt.zmax-0.02)then
         !Strat corenum
         conc(i, 1) = 1.d-6 !nad
         conc(i, 2) = 1.d-4 !nadh
         conc(i, 3) = 0.d0 !fad
         conc(i, 4) = 0.d0 !riboflavin
         conc(i, 5) = 0.d0 !tyrosine
         conc(i, 6) = 0.d0 !trytophan
      elseif(z.gt.zmax-0.1)then
         !Living Epidermis
         conc(i, 1) = 1.d-3 !nad
         conc(i, 2) = 1.d-4 !nadh
         conc(i, 3) = 0.d0 !fad
         conc(i, 4) = 0.d0 !riboflavin
         conc(i, 5) = 1.d-3 !tyrosine
         conc(i, 6) = 0.d0 !trytophan 
      elseif(z.gt.zmax-0.28)then
         !PaPillary Dermis
         conc(i, 1) = 0.d0 !nad
         conc(i, 2) = 0.d0 !nadh
         conc(i, 3) = 0.d0 !fad
         conc(i, 4) = 5.d-4 !riboflavin
         conc(i, 5) = 0.d0 !tyrosine
         conc(i, 6) = 1.d-5 !trytophan
      elseif(z.gt.zmax-2.1)then
         !Reticular Dermis    
         conc(i, 1) = 0.d0 !nad
         conc(i, 2) = 0.d0 !nadh
         conc(i, 3) = 0.d0 !fad
         conc(i, 4) = 0.d0 !riboflavin
         conc(i, 5) = 0.d0 !tyrosine
         conc(i, 6) = 0.d0 !trytophan     
      elseif(z.lt.zmax-2.1)then
         !Hypodermis
         conc(i, 1) = 0.d0 !nad
         conc(i, 2) = 0.d0 !nadh
         conc(i, 3) = 0.d0 !fad
         conc(i, 4) = 0.d0 !riboflavin
         conc(i, 5) = 0.d0 !tyrosine
         conc(i, 6) = 0.d0 !trytophan
      end if
   end do
   call opt_set()

   !****************** Calculate equatorial and polar optical depths ****
   taueq1=0.
   taupole1=0.
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
