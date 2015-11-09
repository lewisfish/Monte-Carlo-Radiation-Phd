      subroutine reader(hgg,mua,mus,opt_parmas,cnt)
      
      implicit none
      
      character(len=100) :: opt_parmas,line
      double precision,intent(out) ::hgg(8),mua(8),mus(8)
      integer :: io,cnt,i
      
      open(1,file=opt_parmas,status='old',iostat=io)
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
      open(2,file=opt_parmas,status='old',iostat=io)

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
      
      character(*),intent(in) :: line
      character(len=30) :: cfg_param
      integer :: cnt,i
      double precision :: hgg(cnt),mua(cnt),mus(cnt)
      
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
      end subroutine reader
