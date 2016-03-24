subroutine findsmax(xmax,ymax,zmax,xcur,ycur,zcur,delta,smax,tflag)

use photon_vars, only : nxp,nyp,nzp

implicit none

DOUBLE PRECISION xcur,ycur,zcur,dsx,dsy,dsz,delta,smax
DOUBLE PRECISION xmax,ymax,zmax
logical tflag

dsx=0.
dsy=0.
dsz=0.

if(nxp.gt.0.) then
   dsx=(2.*xmax-xcur)/nxp
elseif(nxp.lt.0.) then
   dsx=-xcur/nxp
elseif(nxp.eq.0.) then
   dsx=1.e2*xmax
endif

if(nyp.gt.0.) then
   dsy=(2.*ymax-ycur)/nyp
elseif(nyp.lt.0.) then
   dsy=-ycur/nyp
elseif(nyp.eq.0.) then
   dsy=1.e2*ymax
endif

if(nzp.gt.0.) then
   dsz=(2.*zmax-zcur)/nzp
elseif(nzp.lt.0.) then
   dsz=-zcur/nzp
elseif(nzp.eq.0.) then
   dsz=1.e2*zmax
endif

smax=min(dsx,dsy,dsz)
if(smax.lt.delta) then
   tflag=.TRUE.
   return
endif
end subroutine findsmax
