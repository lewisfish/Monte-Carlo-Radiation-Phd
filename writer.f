      subroutine writer(image,reflc,trans,Nbins,fileplace,flu,deposit,
     +                  cbinsnum,fluro,cnt)
      
      implicit none
      
      integer i,j,Nbins,cbinsnum,cnt
      real reflc(1:cnt,1:cnt),trans(1:cnt,1:cnt)
      real image(-((Nbins-1)/2):((Nbins-1)/2),
     + -((Nbins-1)/2):((Nbins-1)/2))
      real flu(-((Nbins-1)/2):((Nbins-1)/2),-((
     +      Nbins-1)/2):((Nbins-1)/2))
      real deposit(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum)
      real fluro(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum)
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

      open(68,file=fileplace//'deposit.dat')
      open(67,file=fileplace//'fludeposit.dat')
      do i=0,cbinsnum
      
            write(68,*) (deposit(i,100,j),j=0,cbinsnum)
            write(67,*) (fluro(i,100,j),j=0,cbinsnum)
      
      end do
      
      open(69,file=fileplace//'image.dat')
      open(70,file=fileplace//'flu.dat')
      do i=-(Nbins-1)/2,(Nbins-1)/2 - 1
      
            write(69,*) (image(j,i),j=-(Nbins-1)/2,(Nbins-1)/2 -  1)
            write(70,*) (flu(i,j),j=-(Nbins-1)/2,(Nbins-1)/2 - 1)
      end do
      
      close(67)
      close(68)
      close(69)
      close(70)
      
      return
      end
