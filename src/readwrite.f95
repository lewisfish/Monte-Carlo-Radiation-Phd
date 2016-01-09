program rea
implicit none

c**** Program reads in unformatted 3D arrays output by 
c**** my grid code. Outputs 2D slices through the grids
c**** which can then be read and displayed by gnuplot, 
c**** Mathematica, etc

integer i,j,Nbins,nxg,cbinsnum
parameter(nxg=201,cbinsnum=200,Nbins=401)
real*8 jmean(nxg+3,nxg+3,nxg+3,4)
real*8 deposit(cbinsnum,cbinsnum,cbinsnum,4)
real*8 image(-((Nbins-1)/2):((Nbins-1)/2),
     +      -((Nbins-1)/2):((Nbins-1)/2),4)
character(*),parameter::
     +      fileplace="/home/st-andrews/lm959/data/"


c**** Read in jmean grid as unformatted array
open(10,file=fileplace//'jmean/jmean.dat',status='unknown',
     +        form='unformatted',access="stream")
     read(10) jmean
close(10)
c**** Read in deposit grid as unformatted array
open(11,file=fileplace//'deposit/deposit.dat',status='unknown',
     +        form='unformatted')
     read(11) deposit
close(11)

c**** Read in image grid as unformatted array
open(9,file=fileplace//'im/images.dat',status='unknown',
     +        form='unformatted')
     read(9) image
close(9)

c**** Write out a 2D slice for jmean.     
open(12,file=fileplace//'jmean/809.dat',status='unknown')      
open(13,file=fileplace//'jmean/1064.dat',status='unknown')      
open(14,file=fileplace//'jmean/900.dat',status='unknown')      
open(15,file=fileplace//'jmean/1300.dat',status='unknown')      
do i=1,nxg
write(12,*) (jmean(i,100,j,1), j=1,nxg)
write(13,*) (jmean(i,100,j,2), j=1,nxg)
write(14,*) (jmean(i,100,j,3), j=1,nxg)
write(15,*) (jmean(i,100,j,4), j=1,nxg)
end do
close(12)
close(13)
close(14)
close(15)

c**** Write out a 2D slice for deposit.     

open(16,file=fileplace//'deposit/809.dat',status='unknown')      
open(17,file=fileplace//'deposit/1064.dat',status='unknown')      
open(18,file=fileplace//'deposit/900.dat',status='unknown')      
open(19,file=fileplace//'deposit/1300.dat',status='unknown')      
do i=1,cbinsnum
write(16,*) (deposit(i,100,j,1), j=1,cbinsnum)
write(17,*) (deposit(i,100,j,2), j=1,cbinsnum)
write(18,*) (deposit(i,100,j,3), j=1,cbinsnum)
write(19,*) (deposit(i,100,j,4), j=1,cbinsnum)
end do
close(12)
close(13)
close(14)
close(15)

c**** Write out images. 
    
open(20,file=fileplace//'im/809.dat',status='unknown')
open(21,file=fileplace//'im/1064.dat',status='unknown')      
open(22,file=fileplace//'im/900.dat',status='unknown')      
open(23,file=fileplace//'im/1300.dat',status='unknown')      
do i=-(Nbins-1)/2,(Nbins-1)/2-1
write(20,*) (image(j,i,1),j=-(Nbins-1)/2,(Nbins-1)/2-1)
write(21,*) (image(j,i,2),j=-(Nbins-1)/2,(Nbins-1)/2-1)
write(22,*) (image(j,i,3),j=-(Nbins-1)/2,(Nbins-1)/2-1)
write(23,*) (image(j,i,4),j=-(Nbins-1)/2,(Nbins-1)/2-1)
end do

end
