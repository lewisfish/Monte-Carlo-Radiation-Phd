      subroutine writer(imageGLOBAL,reflc,trans,Nbins,fileplace,
     +            fluGLOBAL,depositGLOBAL,cbinsnum,fluroglobal
     +            ,cnt,jmeanGLOBAL)
      
      implicit none
      
      include 'grid.txt'
      
      integer i,j,Nbins,cbinsnum,cnt
      real reflc(1:cnt,1:cnt),trans(1:cnt,1:cnt)
      real imageGLOBAL(-((Nbins-1)/2):((Nbins-1)/2),
     + -((Nbins-1)/2):((Nbins-1)/2))
      real fluGLOBAL(-((Nbins-1)/2):((Nbins-1)/2),-((
     +      Nbins-1)/2):((Nbins-1)/2))
      real depositGLOBAL(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum)
      real fluroGLOBAL(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum)
      character(*) fileplace 

     
      
      open(61,file=fileplace//'testreflc.dat')
      do i=1,cnt
            write(61,*) (reflc(i,j),j=1,cnt)
      end do
      close(61)
      
      open(62,file=fileplace//'testtrans.dat')
      
      do i=1,cnt
            write(62,*) (trans(i,j),j=1,cnt)
      end do
      close(62)
      
      open(63,file=fileplace//'jmeanflu.dat')
      open(64,file=fileplace//'jmean.dat')
      do i=1,nxg
            write(63,*) (jmeanGLOBAL(i,100,j,2),j=1,nzg)
            write(64,*) (jmeanGLOBAL(i,100,j,1),j=1,nzg)
      end do
      close(63)
      close(64)
      

      open(68,file=fileplace//'deposit.dat')
      open(67,file=fileplace//'fludeposit.dat')
      do i=0,cbinsnum
      
            write(68,*) (depositGLOBAL(i,100,j),j=0,cbinsnum)
            write(67,*) (fluroGLOBAL(i,100,j),j=0,cbinsnum)
      
      end do
      

      open(69,file=fileplace//'image.dat')
      open(70,file=fileplace//'flu.dat')
      do i=-(Nbins-1)/2,(Nbins-1)/2-1
      
            write(69,*) (imageGLOBAL(j,i),j=-(Nbins-1)/2,(Nbins-1)/2-1)
            write(70,*) (fluGLOBAL(j,i),j=-(Nbins-1)/2,(Nbins-1)/2-1)
      end do
      
      open(71,file=fileplace//'rawimage.dat')
      open(72,file=fileplace//'rawfluimage.dat')
      write(71,*) imageGLOBAL
      write(72,*) fluGLOBAL
      
      close(67)
      close(68)
      close(69)
      close(70)
      close(71)
      close(72)
      
      return
      end
