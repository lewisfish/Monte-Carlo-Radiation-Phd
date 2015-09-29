      subroutine writer(cnt,reflc,trans,image,Nbins,fileplace)
      
      implicit none
      
      integer i,j,cnt,Nbins
      real reflc(1:cnt,1:cnt),trans(1:cnt,1:cnt)
      real image(-((Nbins-1)/2):((Nbins-1)/2),
     + -((Nbins-1)/2):((Nbins-1)/2))
           character(*) fileplace 

     
      
      open(67,file=fileplace//'testreflc.dat')
      do i=1,cnt
            write(67,*) (reflc(i,j),j=1,cnt)
      end do
      close(67)
      
      open(68,file=fileplace//'testtrans.dat')
      
      do i=1,cnt
            write(68,*) (trans(i,j),j=1,cnt)
      end do
      close(68)
      
      open(69,file=fileplace//'image.dat')
      
      do i=-(Nbins-1)/2,(Nbins-1)/2
            write(69,*) (image(j,i),j=
     +      -(Nbins-1)/2,(Nbins-1)/2)
      end do
      
      return
      end
