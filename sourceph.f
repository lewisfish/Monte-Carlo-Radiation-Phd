      subroutine sourceph(xp,yp,zp,nxp,nyp,nzp,
     +                    sint,cost,sinp,cosp,phi,
     +                    xmax,ymax,zmax,twopi,
     +                    xcell,ycell,zcell,nxg,nyg,nzg,iseed)

      implicit none

      include 'photon.txt'

      integer xcell,ycell,zcell,nxg,nyg,nzg,iseed
      real xmax,ymax,zmax,w,twopi
      real ran2

c***** emit photon isotropically from origin
!      xp=0.
!      yp=0.
!      zp=zmax
!      w=0.2

!      do while(xp**2+yp**2 .gt. w**2)
!            xp=2.*ran2(iseed)-1.
!            yp=2.*ran2(iseed)-1.
!      end do
      
!      cost=2.*ran2(iseed)-1.
!      sint=(1.-cost*cost)
!      if(sint.le.0.)then
!        sint=0.
!      else
!        sint=sqrt(sint)
!      endif
      
!      phi=twopi*ran2(iseed)
!      cosp=cos(phi)
!      sinp=sin(phi)


c**** emit uniformly across surface

      zp=zmax
      xp=2.*xmax*ran2(iseed)-xmax
      yp=2.*ymax*ran2(iseed)-ymax

      phi=0.
      cosp=0.
      sinp=0.
      sint=0.
      cost=-1.

c***** Set photon direction cosines for direction of travel *********
      nxp=sint*cosp  
      nyp=sint*sinp
      nzp=cost

c*************** Linear Grid *************************
      xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
      ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
      zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
c*****************************************************

      return
      end

