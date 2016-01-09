subroutine peelingoff(xmax,ymax,zmax,xcell,ycell,zcell,delta &
                      ,image,Nbins,v,g2,hgg,sintim,costim,sinpim,cospim) 

use constants,only : xface,yface,zface,rhokap,pi
use photon,only : xp,yp,zp,nxp,nyp,nzp

implicit none

real xim,yim,tau1,prob,v(3),sinpim
real zmax,bin_width,hgfact,sintim,xmax
real delta,ymax,cosa,costim,cospim
integer binx,biny,xcell,ycell,zcell,Nbins
real image(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2):((Nbins-1)/2),4),g2(1),hgg(1)

! set bin width for images
bin_width=4.*xmax/Nbins

! calculate distance to grid edge
call taufind1(xmax,ymax,zmax,v,tau1,xcell,ycell,zcell,delta)

!calculate angle the peeled off photon will go in
cosa=nxp*v(1)+nyp*v(2)+nzp*v(3)
 
! calc image pixel peeled off photon will 'hit'     
xim=yp*cospim-xp*sinpim
yim=zp*sintim-yp*costim*sinpim-xp*costim*cospim

!calc bin integers for image
binx=floor(xim/bin_width)
biny=nint(yim/bin_width)



!calc weighting for this direction to happen
!      hgfact=(1.-g2(1))/(4.*pi*(1.+g2(1)-2.*hgg(1)*cosa)**(1.5))
 
!calc total weight of peeled off photon
prob=1.*exp(-tau1)

image(binx,biny,1)=image(binx,biny,1)+prob


end subroutine peelingoff
