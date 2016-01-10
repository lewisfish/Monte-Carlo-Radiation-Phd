subroutine noisey(xcell,ycell,cnt)

use photon, only : nxp,nyp,nzp,cost,sint,cosp,sinp
use iarray, only : noise

implicit none

integer cnt,xcell,ycell
real newNorm(3),Norm(3),xvec(3),yvec(3)
real theta,ygrad,xgrad


!      sets vectors in stupid fashion
xvec(1)=1.
xvec(2)=0.
xvec(3)=0.

yvec(1)=0.
yvec(2)=1.
yvec(3)=0.

norm(1)=0.
norm(2)=0.
norm(3)=1.

!     calculates x and y gradients for bump using finite difference method  
if(xcell.le.1)xcell=2
if(ycell.le.1)ycell=2

xgrad=noise(xcell-1,ycell)+noise(xcell+1,ycell)

ygrad=noise(xcell,ycell-1)+noise(xcell,ycell+1)

if(xgrad.lt.1.E-8)then
newNorm=norm + ygrad*yvec
elseif(ygrad.lt.1.E-8)then
newNorm=norm + xgrad*xvec
elseif((ygrad.lt.1.E-8).and.(xgrad.lt.1.E-8))then
newNorm=norm
else     
!     calculates new norm for bump pixel
newNorm=norm + xgrad*xvec + ygrad*yvec
end if



!    calculate proper theta for fresnel reflec
theta=newNorm(3)/sqrt(newNorm(1)**2+newNorm(2)**2+newNorm(3)**2)


cost=-theta      
sint=sqrt(1.-cost**2)
sinp=0.
cosp=0.

nxp=sint*cosp  
nyp=sint*sinp
nzp=cost
!            
end
