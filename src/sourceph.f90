MODULE sourceph_mod

implicit none
save

CONTAINS
   subroutine sourceph(xcell,ycell,zcell,iseed)

   use constants, only : nxg,nyg,nzg,TWOPI,xmax,ymax,zmax
   use photon_vars

   implicit none


   integer xcell,ycell,zcell,iseed
   DOUBLE PRECISION w,r1,phigauss,zf
   DOUBLE PRECISION xf,temp
   real ran2

   !***** emit photon from a circle on surface
   ! min+ran2*(max-min)
!   zp=zmax
!!!   !radius in mm
!!!   w=0.3
!         xp=0.
!         yp=0.
!   xp=xmax
!   yp=ymax
!!   do while(xp**2+yp**2 .gt. w**2.)
!       xp=w*2.*ran2(iseed)-w
!       yp=w*2.*ran2(iseed)-w
!   end do
   !***** emit isotropically from a point

!         xp=0.
!         yp=0.

   !      cost=2.*ran2(iseed)-1.
   !      sint=(1.-cost*cost)
   !      if(sint.le.0.)then
   !        sint=0.
   !      else
   !        sint=sqrt(sint)
   !      endif

   !      phi=TWOPI*ran2(iseed)
   !      cosp=cos(phi)
   !      sinp=sin(phi)


   !**** emit uniformly across surface

         zp=zmax
         xp=2.*xmax*ran2(iseed)-xmax
         yp=2.*ymax*ran2(iseed)-ymax

   !**** Collimated Gaussian Beam

!         r1=w*sqrt(-log(ran2(iseed)))
!         phigauss=TWOPI*ran2(iseed) 
!         xp=r1*cos(phigauss)
!         yp=r1*sin(phigauss)

   !**** Focused Gaussian Beam
   !      zf=.05
   !      r1=w*sqrt(-log(ran2(iseed)))      
   !      phigauss=TWOPI*ran2(iseed) 
   !      xp=r1*cos(phigauss)
   !      yp=r1*sin(phigauss)
   !      
   !      if(2.*ran2(iseed)-1..lt.0.)then
   !            xf=-1.*w*sqrt(-log(ran2(iseed)))
   !      else
   !            xf=w*sqrt(-log(ran2(iseed)))
   !      end if  
   !      temp=sqrt((xp-xf)**2+zf**2)
   !      sint=(xp-xf)/temp
   !      cost=-zf/temp
   !      cosp=1.
   !      sinp=0.



   !***** Set photon direction cosines for direction of travel *********

   phi=0.
   cosp=0.
   sinp=0.
   sint=0.
   cost=-1.

   nxp=sint*cosp  
   nyp=sint*sinp
   nzp=cost

   !*************** Linear Grid *************************
   xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
   ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
   zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
   !*****************************************************

   end subroutine sourceph
end MODULE sourceph_mod
