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
   integer                       :: id,error,numproc,io
   integer status(MPI_STATUS_SIZE) 

! init mpi
call MPI_init(error)
! get number of processes
call MPI_Comm_size(MPI_COMM_WORLD,numproc,error)
! get individual process id
call MPI_Comm_rank(MPI_COMM_WORLD,id,error)
   
   if(id.eq.0)print*,'Starting N-Models on',numproc,'cores'
   filen = '/home/lewis/Desktop/Ga-simple-monte/N-models/data/'

   call directory(id)
   call reader1


   N = 20   !number of grid points => N**3 points to be run
   
   !black magic i will never remember
   !bascically splits up the loop over all points on to 2 processors each with even amount of points to be run on each
   !i.e running double the amount of photons for each point
   !up is the upper limit for the loop for that core and low is the lower bound
   !dosent hold for all N and numproc :(
!   num = int(2.*real(n)/real(numproc) - 1)
!   
!   do i = 0, numproc,2
!      if(id==i)then
!         if(i==0)then
!            low = i*(real(num)+1./2.)+1
!            up = (num+1)
!         else if(i.gt.2)then
!            low = i*int((real(num)+1.)/2.)+1
!            up = (i-((real(i)/2.)-1))*(num+1)
!         else
!            low = i*int((real(num)+1.)/2.)+1
!            up = (num+1)*i
!         end if
!         call MPI_SEND(up, 1, MPI_INTEGER, i+1,1,MPI_COMM_WORLD,error)
!         call MPI_SEND(low, 1, MPI_INTEGER, i+1,1,MPI_COMM_WORLD,error)
!      end if
!   end do
!   
!   do i=1,numproc,2
!      if(id==i)then
!         call MPI_RECV(up, 1, MPI_INTEGER, i-1, 1,MPI_COMM_WORLD,status,error)
!         call MPI_RECV(low, 1, MPI_INTEGER, i-1, 1,MPI_COMM_WORLD,status,error)
!      end if
!   end do
!!   print*,id,low,up
!   call MPI_Barrier(MPI_COMM_WORLD,error)

   
   allocate(conc(3,N,N,N))
   
!   if(mod(id,2)==0)then   
!    l = id + 11
!      open(l,file='test-'//achar(id+48)//'.pat',iostat=io)
!   end if
!   call cpu_time(start)
!   call MPI_Barrier(MPI_COMM_WORLD, error)
!   do i = low, up
!   
!      f = dble(i)*((1d-2-1d-4)/dble(N))-1d-4
!      do j = 1, n
!            g = dble(j)*((1d-3-1d-5)/dble(N))-1d-5
!         do k = 1, n
!               h = dble(k)*((1d-4-1d-6)/dble(N))-1d-6
!               conc(:,i,j,k) = (/ f, g, h/)


!               call MPI_Barrier(MPI_COMM_WORLD, error)
!               call mcpolar(conc(:,i,j,k), i, j, k, error, numproc, id, .false.)
!               call MPI_Barrier(MPI_COMM_WORLD, error)
!               if(mod(id,2)==0)then
!                  write(l,"(3(I3.3,1X),3E15.8)") i,j,k,conc(:,i,j,k)
!                  print*,i,j,k,id
!               end if
!         end do
!      end do
!      if(id.eq.0)then
!         write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
!                " Percent Complete: ", (((numproc/2d0)*dble(i)/dble(N)))*100.0, "%"
!      end if
!      call MPI_Barrier(MPI_COMM_WORLD, error)
!   end do
!   
!   call MPI_Barrier(MPI_COMM_WORLD, error)
!   if(id==0)print*,' '
!   call cpu_time(finish)
!   if(finish-start.ge.60..and.id==0)then
!      print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
!   elseif(id==0)then
!      print*, 'time taken ~',floor(finish-start/60.),'s'
!   end if
!   
!   if(id==0)then
!   
!      call system('zip -j /home/lewis/phdshizz/grid/data/test.zip /home/lewis/phdshizz/grid/data/run_test/*.dat')
!   end if
!   call MPI_Barrier(MPI_COMM_WORLD, error)
   conc(:,1,1,1) = (/1.d-3,1.d-4,1.d-5/)
   call mcpolar(conc(:,1,1,1), 111, 111, 111, error, numproc, id, .false.) 
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
