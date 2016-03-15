      program rea
      implicit none

c**** Program reads in unformatted 3D arrays output by 
c**** Kenny's grid code. Outputs 2D slices through the grids
c**** which can then be read and displayed by gnuplot, 
c**** Mathematica, etc

      integer i,j,Nbins,nxg,nyg,nzg,cbinsnum,k
      parameter(nxg=200,nyg=200,nzg=400,Nbins=401,cbinsnum=200)
      real jmean(nxg+3,nyg+3,nzg+3,4),fluro_pos(nxg,nyg,nzg)
      real deposit(cbinsnum,cbinsnum,cbinsnum),jm(nzg)
      real image(-((Nbins-1)/2):((Nbins-1)/2),-((
     +      Nbins-1)/2):((Nbins-1)/2),4),trans(nxg,nyg)
      character(len=100) fileplace
      fileplace ='/home/lewis/phdshizz/grid/data/'

      jm=0.
c**** Read in jmean grid as unformatted array
      open(10,file='jmean/jmean.dat',form='unformatted',
     +access='direct',recl=((nxg+3)*(nyg+3)*(nzg+3))*4*4)
           read(10,rec=1) jmean
      close(10)
      
!      open(32,file='trans.dat',form='unformatted',access='direct'
!     +,recl=cbinsnum*cbinsnum*4)
!           read(32,rec=1) trans
!      close(32)

!      open(14,file='im/images.dat',form='unformatted'
!     + ,access='direct',recl=(Nbins**2)*4*4)
!           read(14,rec=1) image
!      close(14)
!      
!      open(18,file='deposit/deposit.dat',
!     + form='unformatted',access='direct'
!     + ,recl=(cbinsnum**3)*4*4)
!           read(18,rec=1) deposit
!      close(18)
!      
       open(15,file='fluro_pos.dat',
     + form='unformatted',access='direct'
     + ,recl=((nxg)*(nyg)*(nzg))*4*4)
           read(15,rec=1) fluro_pos
      close(15)
      

c**** Write out a 2D slice. Change the numerical index
c**** to choose the location of the slice in x, y, or z  

!      open(56,file='transt.dat')
!      do i=1,cbinsnum
!            write(56,*) (trans(i,j) ,j=1,cbinsnum)
!      end do
      open(10,file='jmean/jmeanslice4.dat')      
!      open(11,file='jmean/1064.dat')      
!      open(12,file='jmean/900.dat')      
!      open(13,file='jmean/1300.dat')      
      do i=1,nxg
            write(10,*) (jmean(i,100,j,1), j=1,nzg)
!            write(11,*) (jmean(i,100,j,2), j=1,nzg)
!            write(12,*) (jmean(i,100,j,3), j=1,nzg)
!            write(13,*) (jmean(i,100,j,4), j=1,nzg)
      end do
      close(10)
      close(11)
      close(12)
      close(13)
      
      do k=1,nzg
         do j=1,nyg
            do i=1,nxg
               jm(k)=jm(k)+jmean(i,j,k,1)
            end do
         end do
      end do
      !normalise j=(L/#xcells*#ycells*N*Vol)*total pathlengths
!      jmean=jmean/(2.*500000.*(1./nxg)**3.)
!      jm=(jm)/(nxg*nyg*2.*500000.*(1./nxg)**3.)
       jm=jm/(nxg*nyg)
      
      open(56, file='jmean/jmean430.dat')
      open(58, file='jmean/jmeancross.dat')

      do j=nzg,int(nzg/2),-1
         write(56,*) ,real(nzg-j)*(2./nzg),jm(j)
         write(58,*) (jmean(100,i,j,1),i=1,nyg)
      end do
      close(11)
      
!      open(18,file='deposit/809.dat')      
!      open(19,file='deposit/1064.dat')      
!      open(20,file='deposit/900.dat')      
!!      open(21,file='deposit/1300.dat')   
!      do i=1,cbinsnum
!            write(18,*) (deposit(i,j,100), j=1,cbinsnum)
!            write(19,*) (deposit(i,100,j), j=1,cbinsnum)
!            write(20,*) (deposit(100,i,j), j=1,cbinsnum)
!!            write(21,*) (deposit(i,j,k), j=1,cbinsnum)
!      end do
!      close(18)
!      close(19)
!      close(20)
!      close(21)
!      
!      open(14,file='im/809.dat')      
!      open(15,file='im/1064.dat')      
!      open(16,file='im/900.dat')      
!      open(17,file='im/1300.dat')      
!      do i=-((Nbins-1)/2),((Nbins-1)/2)
!      write(14,*) (image(i,j,1), j=-((Nbins-1)/2),((Nbins-1)/2))
!      write(15,*) (image(i,j,2), j=-((Nbins-1)/2),((Nbins-1)/2))
!      write(16,*) (image(i,j,3), j=-((Nbins-1)/2),((Nbins-1)/2))
!      write(17,*) (image(i,j,4), j=-((Nbins-1)/2),((Nbins-1)/2))
!      end do
!      close(14)
!      close(15)
!      close(16)
!      close(17)
!      
      open(19,file='fluro_postion.dat')
      do i=1,nxg
         write(19,*) (fluro_pos(100,i,j), j=1,nzg)
      end do
      close(19)
      
      return
      end
