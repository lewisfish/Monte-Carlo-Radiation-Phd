MODULE writer_mod

implicit none
save

CONTAINS
   subroutine writer(nphotons, numproc)

   use constants,only : nxg,nyg,nzg,fileplace,Nbins,cbinsnum,xmax,ymax,zmax
   use iarray,only : jmeanGLOBAL,fluroexitGLOBAL!,imageGLOBAL,depGLOBAL,depositGLOBAL,transGLOBAL &
                    ! ,fluro_posGLOBAL,followGLOBAL

   implicit none
   
   character(len=4), dimension(7) :: fn
   integer :: j,irec,i,nphotons, numproc
   
   fn=(/'nadh','ribo','try ','tyro','fad ','nad ','tot '/)
   
   jmeanGLOBAL=jmeanGLOBAL * (1.d0/(numproc*nphotons*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))
   
   INQUIRE(iolength = irec) jmeanGLOBAL
   open(62,file=trim(fileplace)//'jmean/jmean 635.dat',access='direct',form='unformatted',recl=irec)
   write(62,rec=1) jmeanGLOBAL
   close(62)

!   open(78,file=trim(fileplace)//'jmean/test.dat')
!   do i=1,nxg
!      write(78,*) (jmeanGLOBAL(i,50,j,1),j=1,nzg)
!   end do
!   open(64,file=trim(fileplace)//'deposit/deposit.dat',access='direct',form='unformatted',recl=(cbinsnum**3)*4*4)
!   write(64,rec=1) depositGLOBAL
!   close(64)

!   open(65,file=trim(fileplace)//'deposit/f-dep.dat')
!   do j=1,cbinsnum
!      write(65,*) depGLOBAL(j)
!   end do
!   close(64)

!   print*,'written out deposit'

!   open(69,file=trim(fileplace)//'im/images.dat',access='direct',form='unformatted',recl=(Nbins**2)*4*4)
!   write(69,rec=1) imageGLOBAL
!   close(69)

!   open(70,file=trim(fileplace)//'trans.dat',access='direct',form='unformatted',recl=cbinsnum*cbinsnum*4)  
!   write(70,rec=1) transGLOBAL
!   close(70)

   do i=1,7
   open(71,file=trim(fileplace)//trim(fn(i))//'fluro.dat')
   do j=1,1000
      write(71,*) fluroexitGLOBAL(j,i)!/real(maxval(fluroexitGLOBAL))
   end do
      close(71)
   end do

   
!   open(72,file=trim(fileplace)//'fluro_pos.dat',access='direct',form='unformatted',recl=((nxg)*(nyg)*(nzg))*4)
!   write(72,rec=1) fluro_posGLOBAL
!   close(72)

   
   end subroutine writer
end MODULE writer_mod
