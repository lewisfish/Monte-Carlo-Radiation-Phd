MODULE subs

implicit none
save

CONTAINS
   SUBROUTINE directory(id)
!  subroutine defines vars to hold paths to various folders   
!   
!   
   use constants,only : cwd,homedir,fileplace,resdir
#ifdef intel
   use ifport
#endif
   implicit none

   integer :: io, id

   !get current working directory
#ifdef intel
   io = getcwd(cwd)
#else
   CALL getcwd(cwd)
#endif
   !get 'home' dir from cwd
   homedir=trim(cwd(1:len(trim(cwd))-3))
   !get data dir
   fileplace=trim(homedir)//'data/'
   
   !checks to see if data folder exists, if not creates it.
#ifdef intel
   io = chdir(fileplace)
#else
   call chdir(fileplace,io)
#endif
   if(io.ne.0.and.id.eq.0)then
      print*,'data directory does not exist...',id
      print*, 'creating directory...'
#ifdef intel
      io = system("mkdir "//fileplace)
      io = chdir(fileplace)
      io = system("mkdir jmean/")
      io = system("mkdir deposit/")
      io = system("mkdir im/")
#else
      call system("mkdir "//fileplace)
      call chdir(fileplace,io)
      call system("mkdir jmean/")
      call system("mkdir deposit/")
      call system("mkdir im/")
#endif
      print*, 'created directory ',trim(fileplace)
   end if
   
   !get res dir
   resdir=trim(homedir)//'res/'
   
   end SUBROUTINE directory
   
   SUBROUTINE zarray
   
   use iarray
   
   !sets all arrays to zero
   implicit none
   
!   reflc=0.
!   trans=0.
!   image=0.
!   deposit=0.
!   dep=0.
   jmean=0.
!   transGLOBAL=0.
!   jmeanGLOBAL=0.
!   imageGLOBAL=0.
!   depositGLOBAL=0.
!   depGLOBAL=0.
   xface=0.
   yface=0.
   zface=0.
   rhokap=0.
   rgb=0.d0
!   rgbGLOBAL = 0.0d0
   conc = 0.d0
   fluroexit=0
!   fluroexitGLOBAL=0
!   fluro_pos=0.
!   fluro_posGLOBAL=0.
   
   end SUBROUTINE zarray

   SUBROUTINE alloc_array(flag)
!  subroutine allocates allocatable arrays
!   
!   
   use iarray
   use constants,only : nxg,nyg,nzg,Nbins,cbinsnum
   
   implicit none
   
   integer :: flag
   
   if(flag .eq. 0)then!allocate local arrays
   
   allocate(xface(nxg+3),yface(nyg+3),zface(nzg+3))
   allocate(rhokap(nxg+3,nyg+3,nzg+3,1),albedo(nxg+3,nyg+3,nzg+3,1))
   allocate(rgb(nxg,nyg,3), conc(nzg, 6))
   
!   allocate(image(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2):((Nbins-1)/2),4))
!   allocate(reflc(cbinsnum,cbinsnum),dep(cbinsnum))
!   allocate(fluro_pos(nxg,nyg,nzg),fluro_posGLOBAL(nxg,nyg,nzg))
   
!   allocate(imageGLOBAL(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2):((Nbins-1)/2),4),depGLOBAL(cbinsnum))
!   allocate(deposit(cbinsnum,cbinsnum,cbinsnum),depositGLOBAL(cbinsnum,cbinsnum,cbinsnum))
!   allocate(transGLOBAL(cbinsnum,cbinsnum),trans(cbinsnum,cbinsnum))
   allocate(jmean(nxg+3,nyg+3,nzg+3,4))
      allocate(fluroexit(1000,7))
   else!deallocate local arrays and allocate GLOBAL
      allocate(jmeanGLOBAL(nxg+3,nyg+3,nzg+3,4), rgbGLOBAL(nxg,nyg,3),fluroexitGLOBAL(1000,7))
         rgbGLOBAL = 0.0d0
         jmeanGLOBAL = 0.0d0
            fluroexitGLOBAL=0
      deallocate(conc,xface,yface,albedo,rhokap)

   end if
   
   end SUBROUTINE alloc_array
   
   
   subroutine sample(array,cdf,wave,iseed)
!      
!  samples a random value from an array based upon its cdf     
!      
      implicit none
      
      integer, intent(IN)                :: iseed
      DOUBLE PRECISION,    intent(IN)    :: array(:, :), cdf(:)
      DOUBLE PRECISION,    intent(OUT)   :: wave
      
      real :: ran2
      DOUBLE PRECISION :: value
      integer :: nlow
      
      value = ran2(iseed)
      
      call search_1D(cdf,nlow,value)
      call lin_inter_1D(array,cdf,value,nlow,wave)
   
   end subroutine sample
   
   subroutine lin_inter_1D(array,cdf,value,nlow,y)
!
!  linear interpolates between values for an array and its cdf
!   
      implicit none
   
      DOUBLE PRECISION,    intent(OUT)  :: y
      DOUBLE PRECISION,    intent(IN)   :: value, array(:, :), cdf(:)
      integer, intent(IN)   :: nlow
   
      y = array(nlow+1, 1) + (array(nlow+2, 1) - array(nlow + 1, 1)) * (value - cdf(nlow))/(cdf(nlow+1) - cdf(nlow))
   
   end subroutine lin_inter_1D
   
   subroutine lin_inter_2D(array,value,length,nlow,y)
!
!  linear interpolation for an array
!
      implicit none

      DOUBLE PRECISION,    intent(OUT)  :: y
      integer, intent(IN)   :: length
      DOUBLE PRECISION,    intent(IN)   :: value, array(length,2)
      integer, intent(IN)   :: nlow
   
      y = array(nlow, 2) + (array(nlow + 1, 2) - array(nlow, 2)) * (value - array(nlow, 1)) / &
                                                         (array(nlow + 1, 1) - array(nlow, 1))
   
   end subroutine lin_inter_2D
   
   subroutine search_1D(array,nlow,value)
!
!  search by bisection for 1D array
!
      implicit none
      
      integer                          :: nup,middle
      integer,             intent(OUT) :: nlow
      DOUBLE PRECISION,    intent(in)  :: array(:), value
      
      nup = size(array)
      nlow = 1
      middle = int((nup + nlow) / 2.d0)

      do while((nup - nlow) .gt. 1)
         middle = int((nup + nlow) / 2.d0)
         if(value .gt. array(middle))then
            nlow = middle
         else
            nup = middle   
         end if
      end do
   end subroutine search_1D
   
   subroutine search_2D(length,array,nlow,value)
!
!  search by bisection for 2D array
!
      implicit none
      
      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      DOUBLE PRECISION,    intent(in)  :: array(length,2),value
      
      nup = length
      nlow = 1
      middle = int((nup+nlow)/2.)

      do while((nup - nlow).gt.1)
         middle = int((nup + nlow)/2.)
         if(value.gt.array(middle,1))then
            nlow = middle
         else
            nup = middle   
         end if
      end do
   end subroutine search_2D
   
   subroutine mk_cdf(array,cdf,length)
!
!  subroutine that creates cdf for an array of values.
!
      implicit none

      integer, intent(IN)    :: length
      DOUBLE PRECISION,    intent(IN)    :: array(length,2)
      DOUBLE PRECISION,    intent(INOUT) :: cdf(length)
      DOUBLE PRECISION                   :: summ
      integer                :: i,j
   
      do j=1,length-1
         summ=0.
         do i=1,j   
            summ=summ+0.5*(array(i+1,2)+array(i,2))*(array(i+1,1)-array(i,1))
         end do
         cdf(j)=summ      
      end do
      cdf=cdf/cdf(length-1)
   
   end subroutine mk_cdf
   
end MODULE subs
