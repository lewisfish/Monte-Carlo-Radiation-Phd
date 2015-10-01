      subroutine peelingoff(xmax,ymax,zmax,cur,
     +                  xface,yface,zface,rhokap,
     +                  xcell,ycell,zcell,delta,xp,yp,zp,
     +                  pi,image,Nbins,flu)
      
      implicit none
      
      real xim,yim,tau1,prob,v(3),phi,theta,pi
      real xp,yp,zp,zmax,bin_width,x1,y1,z1
      real xmax,xface,yface,zface,smax
      real delta,dist,ymax,rhokap,dist2
      integer binx,biny,xcell,ycell,zcell,Nbins,cur
      real image(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2)
     +      :((Nbins-1)/2))
      real flu(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2)
     +      :((Nbins-1)/2))


      bin_width=2.1*xmax/Nbins
      !set image location
      phi=0.
      theta=0.   
            
      x1=xp
      y1=yp
      z1=zp
      
      call gridedge(x1,y1,z1,zmax,dist,v,dist2)

      call taufind1(xp,yp,zp,xmax,ymax,zmax,v,
     +                    xface,yface,zface,rhokap,smax,
     +                    tau1,xcell,ycell,zcell,delta,cur)
      prob=exp(-tau1)/(4.*pi*dist2**2)
      
      xim=y1*cos(phi)-x1*sin(phi)
      yim=z1*sin(theta)-y1*cos(theta)*sin(phi)-x1*cos(theta)*cos(phi)
      
      binx=floor(xim/bin_width)
      biny=nint(yim/bin_width)
      if(cur.eq.1)then
            image(binx,biny)=image(binx,biny)+prob
      else
            flu(binx,biny)=flu(binx,biny)+prob
      end if
      end subroutine peelingoff
