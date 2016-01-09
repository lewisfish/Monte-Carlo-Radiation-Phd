subroutine writer(imageGLOBAL,Nbins,depGLOBAL,depositGLOBAL,cbinsnum,transGLOBAL)

use constants,only : nxg,nyg,nzg,jmeanGLOBAL,fileplace

implicit none

integer Nbins,cbinsnum,j
real imageGLOBAL(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2):((Nbins-1)/2),4)
real transGLOBAL(cbinsnum,cbinsnum)
real depositGLOBAL(cbinsnum,cbinsnum),depGLOBAL(cbinsnum)

open(62,file=fileplace//'jmean/jmean.dat',access='direct',form='unformatted',recl=((nxg+3)*(nyg+3)*(nzg+3))*4*4)

write(62,rec=1) jmeanGLOBAL
close(62)

open(64,file=fileplace//'deposit/deposit.dat',access='direct',form='unformatted',recl=(cbinsnum**2)*4)

write(64,rec=1) depositGLOBAL

close(64)

open(65,file=fileplace//'deposit/f-dep.dat')
do j=1,cbinsnum
write(65,*) depGLOBAL(j)
end do
close(64)

print*,'written out deposit'

open(69,file=fileplace//'im/images.dat',access='direct',form='unformatted',recl=(Nbins**2)*4*4)

write(69,rec=1) imageGLOBAL
close(69)

open(70,file=fileplace//'trans.dat',access='direct',form='unformatted',recl=cbinsnum*cbinsnum*4)
     
write(70,rec=1) transGLOBAL

close(70)
       
return
end
