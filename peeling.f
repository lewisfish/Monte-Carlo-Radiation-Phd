      subroutine peelingoff(xmax,ymax,zmax,cur,nxp,nyp,nzp
     +                  ,xface,yface,zface,rhokap,
     +                  xcell,ycell,zcell,delta,xp,yp,zp,
     +                  pi,image,Nbins,v,g2,hgg,
     +                  sintim,costim,sinpim,cospim,tauflag)
      
            
      implicit none
      
      double precision xim,yim,tau1,prob,v(3),pi,sinpim
      double precision xp,yp,zp,zmax,bin_width,hgfact,sintim
      double precision xmax,xface,yface,zface,smax,nxp,nyp,nzp
      double precision delta,ymax,rhokap,cosa,costim,cospim
      integer binx,biny,xcell,ycell,zcell,Nbins,cur
      logical tauflag
      double precision image(-((Nbins-1)/2):((Nbins-1)/2),
     +      -((Nbins-1)/2):((Nbins-1)/2),4),g2(8),hgg(8)

      ! set bin width for images
      bin_width=4.*xmax/Nbins
      
      ! calculate distance to grid edge
      call taufind1(xp,yp,zp,xmax,ymax,zmax,v,
     +                    xface,yface,zface,rhokap,smax,
     +                    tau1,xcell,ycell,zcell,delta,cur)
     
      !calculate angle the peeled off photon will go in
      cosa=nxp*v(1)+nyp*v(2)+nzp*v(3)
 
      ! calc image pixel peeled off photon will 'hit'     
      xim=yp*cospim-xp*sinpim
      yim=zp*sintim-yp*costim*sinpim-xp*costim*cospim

      !calc bin integers for image
      binx=floor(xim/bin_width)
      biny=nint(yim/bin_width)
      
      ! bin photon in appro image, either fluro or diffuse
      if(cur.eq.1)then
      
            !calc weighting for this direction to happen
      hgfact=(1.-g2(cur))/(4.*pi*(1.+g2(cur)-2.*hgg(cur)*cosa)**(1.5))
 
            !calc total weight of peeled off photon
            prob=hgfact*exp(-tau1)
            
            image(binx,biny,cur)=image(binx,biny,cur)+prob
      else
      
            if(tauflag.eqv..FALSE.)then
      hgfact=(1.-g2(cur))/(4.*pi*(1.+g2(cur)-2.*hgg(cur)*cosa)**(1.5))
            else
                  hgfact=1./(4.*pi)
                  tauflag=.FALSE.
            end if

            !calc total weight of peeled off photon
            prob=hgfact*exp(-tau1)
            image(binx,biny,cur)=image(binx,biny,cur)+prob
      end if
      end subroutine peelingoff
