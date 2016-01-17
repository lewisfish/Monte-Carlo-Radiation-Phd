MODULE subs

implicit none
save

CONTAINS
   SUBROUTINE directory
!  subroutine defines vars to hold paths to various folders   
!   
!   
   use constants,only : cwd,homedir,fileplace,resdir
   implicit none

   integer :: io

   !get current working directory
   CALL getcwd(cwd)

   !get 'home' dir from cwd
   homedir=trim(cwd(1:index(trim(cwd),'/bin')))
   !get data dir
   fileplace=trim(homedir)//'data/'
   
   !checks to see if data folder exists, if not creates it.
   call chdir(fileplace,io)
   if(io.ne.0)then
      print*,'directory does not exist'
      print*, 'creating directory'
      call system("mkdir"//fileplace)
   end if
   
   !get res dir
   resdir=trim(homedir)//'res/'
   
   end SUBROUTINE directory
   
   SUBROUTINE zarray
   
   use iarray
   
   !sets all arrays to zero
   implicit none
   

   reflc=0.
   trans=0.
   image=0.
   deposit=0.
   dep=0.
   jmean=0.
   transGLOBAL=0.
   jmeanGLOBAL=0.
   imageGLOBAL=0.
   depositGLOBAL=0.
   depGLOBAL=0.
   xface=0.
   yface=0.
   zface=0.
   rhokap=0.
   
   end SUBROUTINE zarray

   SUBROUTINE alloc_array
!  subroutine allocates allocatable arrays
!   
!   
   use iarray
   use constants,only : nxg,nyg,nzg,Nbins,cbinsnum
   
   implicit none
   
   allocate(xface(nxg+3),yface(nyg+3),zface(nzg+3))
   allocate(rhokap(nxg+3,nyg+3,nzg+3,4))
   
   allocate(image(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2):((Nbins-1)/2),4))
   allocate(reflc(cbinsnum,cbinsnum),dep(cbinsnum))
   
   allocate(imageGLOBAL(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2):((Nbins-1)/2),4),depGLOBAL(cbinsnum))
   allocate(deposit(cbinsnum,cbinsnum),depositGLOBAL(cbinsnum,cbinsnum))
   allocate(transGLOBAL(cbinsnum,cbinsnum),trans(cbinsnum,cbinsnum))
   allocate(jmean(nxg+3,nyg+3,nzg+3,4)) 
   allocate(jmeanGLOBAL(nxg+3,nyg+3,nzg+3,4))
   
   
   end SUBROUTINE alloc_array
end MODULE subs
