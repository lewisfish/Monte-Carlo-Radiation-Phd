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
   
   subroutine reader1
   
   use iarray,    only : mua_array,mus_array,f_cdf,fluro_array, &
                         e_cdf,excite_array
   use constants, only : resdir
   
   implicit none
   
   integer :: cnt,io,i,j
   
   open(10,file=trim(resdir)//'mua.dat')
   cnt=0
   do

      read(10,*,IOSTAT=io)
     
      if (io < 0) then
         close(10)
         allocate(mua_array(cnt-1,2))
         mua_array=0.
         exit
      else
         cnt=cnt+1
      end if

   end do
   open(20,file=trim(resdir)//'mua.dat') 

   do i=1,cnt-1
      read(20,*) mua_array(i,1),mua_array(i,2)
   end do
   close(20)
   !mus
   open(30,file=trim(resdir)//'interlipid_scat2.dat')
   cnt=0
   do

      read(30,*,IOSTAT=io)
     
      if (io < 0) then
         close(30)
         allocate(mus_array(1:cnt-1,1:2))
         mus_array=0.
         exit
      else
         cnt=cnt+1
      end if

   end do
   
   open(40,file=trim(resdir)//'interlipid_scat2.dat') 
   do i=1,cnt-1
      read(40,*) mus_array(i,1),mus_array(i,2)
   end do
   close(40)
   !fluro
      open(50,file=trim(resdir)//'fluro_louise.dat')
   cnt=0
   do

      read(50,*,IOSTAT=io)
     
      if (io < 0) then
         close(50)
         allocate(fluro_array(cnt,2))
         fluro_array=0.
         exit
      else
         cnt=cnt+1
      end if

   end do
   
   open(50,file=trim(resdir)//'fluro_louise.dat') 
   do i=1,cnt
      read(50,*) fluro_array(i,1),fluro_array(i,2)
   end do
   close(50)
   allocate(f_cdf(cnt))
   f_cdf=0.
   
   !excite
   
   open(90,file=trim(resdir)//'fluro_signal.dat')
   cnt=0
   do

      read(90,*,IOSTAT=io)
     
      if (io < 0) then
         close(90)
         allocate(excite_array(cnt,2))
         excite_array=0.
         exit
      else
         cnt=cnt+1
      end if

   end do
   
   open(11,file=trim(resdir)//'fluro_signal.dat') 
   do i=1,cnt
      read(11,*) excite_array(i,1),excite_array(i,2)
   end do
   close(11)
   allocate(e_cdf(cnt))
   e_cdf=0.

   end subroutine reader1
end MODULE reader_mod
