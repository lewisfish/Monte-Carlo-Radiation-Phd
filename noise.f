      program noise

      implicit none

      real newNorm(3),Norm(3),xvec(3),yvec(3),xgrad,ygrad
      real, allocatable :: noise(:,:)
      integer io,i,j,cnt

            open(13,file='noisedots.dat')
      do
        read(13,*,IOSTAT=io)
        
      if (io < 0) then
            close(13)
            print *, cnt
           allocate(noise(1:cnt-1,1:cnt-1))
           noise=0.
           exit
      else
       cnt=cnt+1
      end if
      end do
      open(14,file='noisedots.dat') 
      do i=1,cnt
            read(14,*) (noise(i,j),j=1,cnt)
      end do
      
      
            xvec=(1.,0.,0.)
            vec=(0.,1.,0.)
            Norm=(0.,0.,1.)      
      
            xgrad=noise(xcur-1.,ycur)+noise(xcur+1.,ycur)
            ygrad=noise(xcur,ycur-1.)+noise(xcur,ycur+1.)
!    calculate new norm for bumpy pixel
            newNorm=norm + xgrad*xvec + ygrad*yvec
            
!    calculate proper theta for fresnel reflec
            theta=180.-acos(1/sqrt(newNorm(1)**2+newNorm(2)**2+newNorm(3)**2))
            print *, theta
      end
