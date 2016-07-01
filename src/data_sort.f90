program simple_GA

   use mpi
   use subs
   use reader_mod
   use main
   implicit none

   real                          :: ran2
   double precision              :: f, g, h, start, finish
   double precision, allocatable :: conc(:,:,:,:)
   integer                       :: i, j, k, N, low, up, num,l
   character(len=255)            :: filen 
   integer                       :: id,error,numproc,io,array(1800,3)
   integer status(MPI_STATUS_SIZE) 

! init mpi
call MPI_init(error)
! get number of processes
call MPI_Comm_size(MPI_COMM_WORLD,numproc,error)
! get individual process id
call MPI_Comm_rank(MPI_COMM_WORLD,id,error)
   
   if(id.eq.0)print*,'Starting N-Models on',numproc,'cores'

id=0
   call directory(id)
   call reader1

   if(id == 0)then
      open(65,file='rosetta_n.txt',iostat=io)
      do i=1,1800
         read(65,*) array(i,:)
!         print*,array(i,:)
      end do
      close(65)
   end if

N=20

!   allocate(conc(3,N,N,N))
!   

   if(id==0)then   
!    l = id + 11
      open(54,file='bob-'//achar(id+48)//'.pat',iostat=io)
      print*,'sdas'
   end if

   call MPI_Barrier(MPI_COMM_WORLD, error)

   do i=1,1800
      f = dble(array(i,1))*((1d-2-1d-4)/dble(N))-1d-4
      g = dble(array(i,2))*((1d-3-1d-5)/dble(N))-1d-5
      h = dble(array(i,3))*((1d-4-1d-6)/dble(N))-1d-6
      
!      conc(:,i,j,k) = (/ f, g, h/)


   call MPI_Barrier(MPI_COMM_WORLD, error)
   call mcpolar((/ f, g, h/), i, j, k, error, numproc, id, .false.)
   call MPI_Barrier(MPI_COMM_WORLD, error)
   if(id==0)then
      write(54,"(3(I3.3,1X),3E15.8)") array(i,:),(/ f, g, h/)
!      print*,i,j,k,id
   end if
   end do
   call MPI_Barrier(MPI_COMM_WORLD, error)
   
   call MPI_Finalize(error)
   Contains 
           
 subroutine readfile(filename, array)
!
!  read in file
!
   double precision, allocatable, dimension(:), intent(inout) :: array
   character(*),                      intent(in)    :: filename
   integer                                          :: cnt, io, i
   
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
            allocate(array(cnt))
            array=0.
            exit
         else
            cnt = cnt + 1
         end if
      end do
      open(20, file = filename) 
      
      !read in data
      do i = 1, cnt 
         read(20, *) array(i)
      end do
      close(20)
   end if
   end subroutine readfile
end program
