      program reader
      
      implicit none
      
      character(len=100) :: filename,line
      real, allocatable :: hgg(:),mua(:),mus(:)
      integer :: io,cnt,i
      
      call get_command_argument(1,filename)
      
      open(1,file=filename,status='old',iostat=io)
      if(io.ne.0)then
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
      open(2,file=filename,iostat=io)
      cnt=cnt-2
      allocate(hgg(cnt),mua(cnt),mus(cnt))
      hgg=0.
      mua=0.
      mus=0.

      do while(.true.)
            i=1
            read(2,'(A)', end=99) line
            print*, line
            
            call linereader(line,cnt,hgg,mus,mua,i)
            i=i+1
      end do
99    continue
      close(2)

      
      contains
      subroutine linereader(line,cnt,hgg,mua,mus,i)
      
      implicit none
      
      character(*),intent(in) :: line
      character(len=30) :: cfg_param
      integer :: cnt,i
      real :: hgg(cnt),mua(cnt),mus(cnt)
      
      cfg_param = adjustl(line(:index(line,'!')-1))
      read(cfg_param,'(F5.3)') hgg(i),mua(i),mus(i) 

      
      end subroutine linereader
      end program reader
