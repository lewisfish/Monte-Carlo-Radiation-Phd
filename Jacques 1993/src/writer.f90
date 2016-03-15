MODULE writer_mod

implicit none
save

CONTAINS
   subroutine writer(xmax,ymax,zmax,nphotons,numproc)

   use constants,only : nxg,nyg,nzg,fileplace,Nbins,cbinsnum
   use iarray,only : jmeanGLOBAL,imageGLOBAL,depGLOBAL,depositGLOBAL,transGLOBAL &
                     ,fluroexitGLOBAL,fluro_posGLOBAL,depthGLOBAL,depthallGLOBAL

   implicit none

   integer j,nphotons,numproc
   real :: xmax,ymax,zmax
   
   jmeanGLOBAL=jmeanGLOBAL * ((2.*xmax)**2./(numproc*nphotons*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))
   
   open(57,file='testf.dat')
   open(56,file='testff.dat')
   open(58,file='all.dat')
   do j=nzg,1,-1
      write(58,*) real(((nzg-j))*real(2.*zmax/nzg)),depthallGLOBAL(j)
      write(57,*) real(((nzg-j))*real(2.*zmax/nzg)),depthGLOBAL(j)
      write(56,*) real(((nzg-j))*real(2.*zmax/nzg)),depthGLOBAL(j)/depthallGLOBAL(j)
   end do
   close(56)
   close(57)
   
   open(62,file=trim(fileplace)//'jmean/jmean.dat',access='direct',form='unformatted',recl=((nxg+3)*(nyg+3)*(nzg+3))*4*4)
   write(62,rec=1) jmeanGLOBAL
   close(62)

   open(64,file=trim(fileplace)//'deposit/deposit.dat',access='direct',form='unformatted',recl=(cbinsnum**3)*4*4)
   write(64,rec=1) depositGLOBAL
   close(64)

   open(65,file=trim(fileplace)//'deposit/f-dep9.dat')
   do j=cbinsnum,1,-1
      write(65,*) depGLOBAL(j),real(((cbinsnum-j))*(1./200.))
   end do
   close(64)

   print*,'written out deposit'

   open(69,file=trim(fileplace)//'im/images.dat',access='direct',form='unformatted',recl=(Nbins**2)*4*4)
   write(69,rec=1) imageGLOBAL
   close(69)

   open(70,file=trim(fileplace)//'trans.dat',access='direct',form='unformatted',recl=cbinsnum*cbinsnum*4)  
   write(70,rec=1) transGLOBAL
   close(70)

   open(71,file=trim(fileplace)//'resfluro8.dat')
   do j=1,1000
      write(71,*) fluroexitGLOBAL(j)!/real(maxval(fluroexitGLOBAL))
   end do
   close(71)
   
   open(72,file=trim(fileplace)//'fluro_pos.dat',access='direct',form='unformatted',recl=((nxg)*(nyg)*(nzg))*4)
   write(72,rec=1) fluro_posGLOBAL
   close(72)

   
   end subroutine writer
end MODULE writer_mod
