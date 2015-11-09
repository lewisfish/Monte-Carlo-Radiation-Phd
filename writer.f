      subroutine writer(imageGLOBAL,reflc,trans,Nbins,fileplace,
     +            fluGLOBAL,depositGLOBAL,cbinsnum,fluroglobal
     +            ,cnt,jmeanGLOBAL)
      
      implicit none
      
      include 'grid.txt'
      
      integer i,j,Nbins,cbinsnum,cnt
      double precision reflc(1:cnt,1:cnt),trans(1:cnt,1:cnt)
      double precision imageGLOBAL(-((Nbins-1)/2):((Nbins-1)/2),
     + -((Nbins-1)/2):((Nbins-1)/2))
      double precision fluGLOBAL(-((Nbins-1)/2):((Nbins-1)/2),
     +      -((Nbins-1)/2):((Nbins-1)/2),4)
      double precision depositGLOBAL(-1:cbinsnum,-1:cbinsnum,-1:
     +cbinsnum)
      double precision fluroGLOBAL(-1:cbinsnum,-1:cbinsnum,-1:
     +cbinsnum,4)
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
      
      open(63,file=fileplace//'jmean/809t1.dat')
      open(64,file=fileplace//'jmean/1064t1.dat')
      open(65,file=fileplace//'jmean/900t1.dat')
      open(66,file=fileplace//'jmean/1300t1.dat')
      do i=1,nxg
            write(63,*) (jmeanGLOBAL(i,100,j,1),j=1,nzg)
            write(64,*) (jmeanGLOBAL(i,100,j,2),j=1,nzg)
            write(65,*) (jmeanGLOBAL(i,100,j,3),j=1,nzg)
            write(66,*) (jmeanGLOBAL(i,100,j,4),j=1,nzg)
      end do
      close(63)
      close(64)
      close(65)
      close(66)
      

      open(67,file=fileplace//'deposit/809t1.dat')
      open(68,file=fileplace//'deposit/1064t1.dat')
      open(73,file=fileplace//'deposit/900t1.dat')
      open(74,file=fileplace//'deposit/1300t1.dat')
      do i=0,cbinsnum
      
            write(67,*) (depositGLOBAL(i,101,j),j=0,cbinsnum)
            write(68,*) (fluroGLOBAL(i,100,j,2),j=0,cbinsnum)
            write(73,*) (fluroGLOBAL(i,100,j,3),j=0,cbinsnum)
            write(74,*) (fluroGLOBAL(i,101,j,4),j=0,cbinsnum)
      
      end do
      

      open(69,file=fileplace//'im/809nmt1.dat')
      open(70,file=fileplace//'im/1064nmt1.dat')
      open(71,file=fileplace//'im/900nmt1.dat')
      open(72,file=fileplace//'im/1300nmt1.dat')      
      do i=-(Nbins-1)/2,(Nbins-1)/2-1
            write(69,*) (imageGLOBAL(j,i),j=-(Nbins-1)/2,(Nbins-1)/2-1)
            write(70,*) (fluGLOBAL(j,i,2),j=-(Nbins-1)/2,(Nbins-1)/2-1)
            write(71,*) (fluGLOBAL(j,i,3),j=-(Nbins-1)/2,(Nbins-1)/2-1)
            write(72,*) (fluGLOBAL(j,i,4),j=-(Nbins-1)/2,(Nbins-1)/2-1)
      end do
       
      close(67)
      close(68)
      close(69)
      close(70)      
      close(71)      
      close(72) 
      close(73)
      close(74)                 
      return
      end
