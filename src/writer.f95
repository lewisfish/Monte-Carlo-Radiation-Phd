MODULE writer_mod

implicit none
save

CONTAINS
   subroutine writer

   use constants,only : nxg,nyg,nzg,fileplace,Nbins,cbinsnum
   use iarray,only : jmean,jmeanGLOBAL,imageGLOBAL,depGLOBAL,depositGLOBAL,transGLOBAL &
                     ,fluroexitGLOBAL

   implicit none

   integer j

   open(62,file=trim(fileplace)//'jmean/jmean.dat',access='direct',form='unformatted',recl=((nxg+3)*(nyg+3)*(nzg+3))*4*4)
   write(62,rec=1) jmeanGLOBAL
   close(62)

   open(64,file=trim(fileplace)//'deposit/deposit.dat',access='direct',form='unformatted',recl=(cbinsnum**2)*4)
   write(64,rec=1) depositGLOBAL
   close(64)

   open(65,file=trim(fileplace)//'deposit/f-dep.dat')
   do j=1,cbinsnum
      write(65,*) depGLOBAL(j)
   end do
   close(64)

   print*,'written out deposit'

   open(69,file=trim(fileplace)//'im/images.dat',access='direct',form='unformatted',recl=(Nbins**2)*4*4)
   write(69,rec=1) imageGLOBAL
   close(69)

   open(70,file=trim(fileplace)//'trans.dat',access='direct',form='unformatted',recl=cbinsnum*cbinsnum*4)  
   write(70,rec=1) transGLOBAL
   close(70)
   
   open(71,file=trim(fileplace)//'resfluro.dat')
   do j=1,1000
      write(71,*) fluroexitGLOBAL(j)
   end do
   close(71)
   
   return
   end
end MODULE writer_mod
