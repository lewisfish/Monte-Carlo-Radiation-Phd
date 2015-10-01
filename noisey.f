      subroutine noisey(xcell,ycell,noise,cnt,angle,nxp,
     +      nyp,nzp,cost,sint,cosp,sinp)

      implicit none
      integer cnt,xcell,ycell
      real newNorm(3),Norm(3),xvec(3),yvec(3),xgrad,ygrad
      real nxp,nyp,nzp,sint,cost,sinp
      real  noise(1:cnt,1:cnt),theta,angle,cosp
      

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

      xgrad=noise(int(xcell-1),int(ycell))+noise(int(xcell+1),
     +      int(ycell))

      ygrad=noise(int(xcell),int(ycell-1))+noise(int(xcell),
     +      int(ycell+1))
!      print *, xgrad,ygrad,xcell
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

!      angle=angle+90.-(180./3.14)*acos(theta)

      cost=-theta      
      sint=sqrt(1.-cost**2)

      nxp=sint*cosp  
      nyp=sint*sinp
      nzp=cost
!            
      end
