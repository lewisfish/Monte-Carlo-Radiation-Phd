      subroutine noisey(xcur,ycur,noise)

      implicit none

      real newNorm(3),Norm(3),xvec(3),yvec(3),xgrad,ygrad,xcur,ycur
      real theta
      real, allocatable :: noise(:,:)
!      integer 
      
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
      
!     calculates x and y gradients for bump usinf finite difference method     
      xgrad=noise(int(xcur-1),int(ycur))+noise(int(xcur+1),int(ycur))
      ygrad=noise(int(xcur),int(ycur-1))+noise(int(xcur),int(ycur+1))
!     calculates new norm for bump pixel
            newNorm=norm + xgrad*xvec + ygrad*yvec
            print *, norm,newNorm
            
!    calculate proper theta for fresnel reflec
      theta=180.-acos(1/sqrt(newNorm(1)**2+newNorm(2)**2+newNorm(3)**2))
      print *, theta
      end
