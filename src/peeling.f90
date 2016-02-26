MODULE peelingoff_mod

implicit none
save

CONTAINS
   subroutine peelingoff(xcell,ycell,zcell,delta &
                         ,v,sintim,costim,sinpim,cospim) 

   use constants, only   : pi,nbins,xmax,ymax,zmax
   use photon_vars, only : xp,yp,zp,nxp,nyp,nzp
   use iarray, only      : xface,yface,zface,rhokap,image
   use opt_prop, only    : hgg,g2
   use taufind
   
   implicit none

   real xim,yim,tau1,prob,v(3),sinpim
   real bin_width,hgfact,sintim
   real delta,cosa,costim,cospim
   integer binx,biny,xcell,ycell,zcell

   ! set bin width for images
   bin_width=4.*xmax/Nbins

   ! calculate distance to grid edge
   call taufind1(v,tau1,xcell,ycell,zcell,delta)

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
end MODULE peelingoff_mod
