      subroutine peelingoff(xmax,ymax,zmax,
     +                  xface,yface,zface,rhokap,d,
     +                  xcell,ycell,zcell,delta,xp,yp,zp,
     +                  sint,sinp,cost,cosp,pi,image,Nbins)
      
      implicit none
      
      real xim,yim,tau1,weight,d,v(3),phi,theta,pi
      real xp,yp,zp,zmax,bin_width,x1,y1,z1
      real xmax,xface,yface,zface
      real delta,sint,sinp,cost,cosp,dist,ymax,rhokap
      integer binx,biny,xcell,ycell,zcell,Nbins
      real image(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2)
     +      :((Nbins-1)/2)),dist2

      bin_width=4.*xmax/Nbins
      !set image location
      phi=0.
      theta=0.   
            
      x1=xp
      y1=yp
      z1=zp
      
      call gridedge(x1,y1,z1,xmax,ymax,zmax,dist,v,dist2)

      call taufind1(xp,yp,zp,xmax,ymax,zmax,v,
     +                    xface,yface,zface,rhokap,
     +                    tau1,xcell,ycell,zcell,delta)

      weight=exp(-tau1)/(4.*pi)
      
      xim=y1*cos(phi)-x1*sin(phi)
      yim=z1*sin(theta)-y1*cos(theta)*sin(phi)-x1*cos(theta)*cos(phi)
      
      binx=floor(xim/bin_width)
      biny=nint(yim/bin_width)

      image(binx,biny)=image(binx,biny)+weight

      end subroutine peelingoff
