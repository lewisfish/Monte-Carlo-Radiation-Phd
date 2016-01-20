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

   open(1,file=trim(resdir)//'opt.params',status='old',iostat=io)
   if(io.ne.0)then
   print*, trim(resdir)//'opt.params','does not exsist'
   print*, 'File cant be opened!'
   print*, 'Exiting...'
   call EXIT(0)
   end if

   cnt=0
   do
   read(1,*,IOSTAT=io)
     
   if (io < 0) then
         close(1)
         exit
   else
         cnt=cnt+1
   end if
   end do
   open(2,file=trim(resdir)//'opt.params',status='old',iostat=io)

   cnt=cnt-1
   !      allocate(hgg(cnt),mua(cnt),mua(cnt))
   hgg=0.
   mua=0.
   mus=0.
   i=1
   do while(.true.)

   read(2,'(A)', end=99) line
   call linereader(line,cnt,hgg,mus,mua,i)
   end do
   99    continue
   close(2)

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
   
   use iarray, only : mua_array,mus_array,fluro_array,cdf
   use constants, only : resdir
   
   implicit none
   
   integer :: cnt,io,i,j
   
   open(1,file=trim(resdir)//'mua.dat')
   cnt=0
   do

      read(1,*,IOSTAT=io)
     
      if (io < 0) then
         close(1)
         allocate(mua_array(cnt-1,2))
         mua_array=0.
         exit
      else
         cnt=cnt+1
      end if

   end do
   open(2,file=trim(resdir)//'mua.dat') 

   do i=1,cnt-1
      read(2,*) mua_array(i,1),mua_array(i,2)
   end do

   !mus
   open(3,file=trim(resdir)//'mus.dat')
   cnt=0
   do

      read(3,*,IOSTAT=io)
     
      if (io < 0) then
         close(3)
         allocate(mus_array(1:cnt-1,1:2))
         mus_array=0.
         exit
      else
         cnt=cnt+1
      end if

   end do
   
   open(4,file=trim(resdir)//'mus.dat') 
   do i=1,cnt-1
      read(4,*) mus_array(i,1),mus_array(i,2)
   end do
   
   !fluro
      open(5,file=trim(resdir)//'fluro.dat')
   cnt=0
   do

      read(5,*,IOSTAT=io)
     
      if (io < 0) then
         close(5)
         allocate(fluro_array(cnt,2))
         fluro_array=0.
         exit
      else
         cnt=cnt+1
      end if

   end do
   
   open(5,file=trim(resdir)//'fluro.dat') 
   do i=1,cnt
      read(5,*) fluro_array(i,1),fluro_array(i,2)
   end do
   
   allocate(cdf(cnt))
   cdf=0.
   open(10,file='/home/lewis/phdshizz/grid/res/testf.dat')
   do i=1,cnt
      write(10,*) fluro_array(i,1),fluro_array(i,2)
   end do
   end subroutine reader1
end MODULE reader_mod
