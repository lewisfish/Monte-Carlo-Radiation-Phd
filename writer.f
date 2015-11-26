      subroutine writer(imageGLOBAL,reflc,trans,Nbins,fileplace,
     +            depositGLOBAL,cbinsnum,cnt,jmeanGLOBAL)
      
      implicit none
      
      include 'grid.txt'
      
      integer i,j,Nbins,cbinsnum,cnt
      double precision reflc(1:cnt,1:cnt),trans(1:cnt,1:cnt)
      double precision imageGLOBAL(-((Nbins-1)/2):((Nbins-1)/2),
     + -((Nbins-1)/2):((Nbins-1)/2),4)
      double precision depositGLOBAL(-1:cbinsnum,-1:cbinsnum,-1:
     +cbinsnum)
      character(*) fileplace
      
      open(63,file=fileplace//'jmean/jmean.dat',form='unformatted')
      write(63) jmeanGLOBAL
      close(63)
      print*,'written Jmean out'

      open(67,file=fileplace//'deposit/deposit.dat',form='unformatted')
      write(67) depositGLOBAL
      close(67)
      print*,'written out deposit'

      

      open(69,file=fileplace//'im/images.dat',form='unformatted')
      write(69) imageGLOBAL

      close(69)
               
      return
      end
