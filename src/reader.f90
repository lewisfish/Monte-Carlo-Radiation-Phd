MODULE reader_mod

implicit none
save

CONTAINS
   subroutine reader0(hgg,mua,mus,cnt)

   use constants,only : resdir

   implicit none

   character(len=100) :: line
   real, intent(out)  :: hgg(1),mua(1),mus(1)
   integer            :: io,cnt,i

   open(10,file=trim(resdir)//'opt.params',status='old',iostat=io)
   if(io.ne.0)then
      print*, trim(resdir)//'opt.params','does not exsist'
      print*, 'File cant be opened!'
      print*, 'Exiting...'
      call EXIT(0)
   end if

   cnt=0
   do
   read(10,*,IOSTAT=io)
     
   if (io < 0) then
         close(10)
         exit
   else
         cnt=cnt+1
   end if
   end do
   open(20,file=trim(resdir)//'opt.params',status='old',iostat=io)

   cnt=cnt-1
   !      allocate(hgg(cnt),mua(cnt),mua(cnt))
   hgg=0.
   mua=0.
   mus=0.
   i=1
   do while(.true.)

      read(20,'(A)', end=99) line
      call linereader(line,cnt,hgg,mus,mua,i)
   end do
   99    continue
   close(20)

   contains
      subroutine linereader(line,cnt,hgg,mua,mus,i)

      implicit none

      character(*), intent(in) :: line
      character(len=30)        :: cfg_param
      integer                  :: cnt,i
      real                     :: hgg(cnt),mua(cnt),mus(cnt)

      if(scan(line,'hgg').eq.0)then
      cfg_param = adjustl(line(:index(line,'!')-1))

      cfg_param = adjustl(cfg_param(:index(cfg_param,' ')+1))
      read(cfg_param,'(F5.3)') hgg(i)
      cfg_param = adjustl(line(index(line,' '):))

      read(cfg_param,'(F5.3)') mus(i)
      cfg_param = adjustl(cfg_param(index(cfg_param,' '):))

      read(cfg_param,'(F5.3)') mua(i)
                  i=i+1
      end if

      end subroutine linereader
   end subroutine reader0
  
!****************************************************************************
   
   subroutine reader1
   
   use iarray,    only : Carotene_array, Bilirubin_array, Oxy_Hb_array, Deoxy_Hb_array, &
                         water_array !, noise
   use constants, only : resdir
   
   implicit none
   
   integer :: cnt

   !Oxy-Hb data
   call readfile_array2D(trim(resdir)//'Oxy-Hb.dat', Oxy_Hb_array, 0, 2)

   !Deoxy-Hb data
   call readfile_array2D(trim(resdir)//'Deoxy-Hb.dat', Deoxy_Hb_array, 0, 2)

   !B-carotene data
   call readfile_array2D(trim(resdir)//'B-carotene.dat', Carotene_array, 0, 2)

   !Bilirubin data
   call readfile_array2D(trim(resdir)//'bilirubin.dat', bilirubin_array, 0, 2)

   !water data
   call readfile_array2D(trim(resdir)//'water absor.dat', water_array, 0, 2)

   !noise data
!   call readfile_array2D(trim(resdir)//'noisedots.dat', noise, 1)

   end subroutine reader1
   
   
!************************************************************************   
  

 subroutine readfile_array2D(filename, array, flag, colsize)
   !
   ! Reads a file to get its length and allocates an array to store the data and reads it in.
   !
   ! subroutine takes filename, the array for dat to be read into, a flag to toggle
   ! between using square array and other shapes, colsize is an optional argument that
   ! specifies the size of the 2nd dimension of the array
   !
   DOUBLE PRECISION, allocatable, dimension(:,:), intent(inout) :: array
   integer,                           intent(in)    :: flag
   integer,                 optional, intent(in)    :: colsize
   character(*),                      intent(in)    :: filename
   integer                                          :: cnt, io, i, j
   
   open(10, file = filename, status = 'OLD', IOSTAT = io)
   if(io .ne. 0)then
      print'(A,A,I2)',filename,' could not be opened. IOSTAT = ',io
      print*,'Exiting...'
      call exit(0)
   else
      cnt = 0
      do       !find file size and allocate array.

         read(10, *, IOSTAT = io)
        
         if (io < 0) then
            close(10)
            if(flag .eq. 0)then
               allocate(array(cnt , colsize))
            elseif(flag .eq. 1)then
               allocate(array(cnt , cnt))
            end if
            array=0.
            exit
         else
            cnt = cnt + 1
         end if

      end do
      open(20, file = filename) 
      
      !read in data
      if(flag .eq. 0)then
         do i = 1, cnt 
            read(20, *) array(i, 1),array(i, 2)
         end do
      elseif(flag .eq. 1)then
            do i = 1, cnt 
            read(20, *) (array(i, j), j = 1, cnt)
         end do
      end if
      close(20)
   end if
   end subroutine readfile_array2D
end MODULE reader_mod
