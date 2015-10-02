      subroutine force(xmax,ymax,zmax,cur,tau,iseed,
     +                xface,yface,zface,rhokap,nxp,nyp,nzp
     +                ,xcell,ycell,zcell,delta,xp,yp,zp)
      
      implicit none
      
      include 'grid.txt'
      
      real xp,yp,zp,smax,tau1,xmax,ymax,zmax,delta
      real v(3),nyp,nzp,nxp,tau,ran2
      integer cur,xcell,ycell,zcell,iseed
      
      v(1)=nxp
      v(2)=nyp
      v(3)=nzp
      
      call taufind1(xp,yp,zp,xmax,ymax,zmax,v,
     +                    xface,yface,zface,rhokap,smax,
     +                    tau1,xcell,ycell,zcell,delta,cur)
     
      tau=-log(1.-ran2(iseed)*(1.-exp(-tau1)))
      
      end subroutine force
