MODULE random

   integer :: iseed

end module

program simple_GA

   use random
   use mpi
   use subs
   use reader_mod
   
   implicit none

   double precision              :: fitval, fitval_m
   real                          :: ran2
   double precision, allocatable :: target_a(:), source(:), m(:), test(:)
   integer                       :: i, j
   
integer          :: error,numproc,id
   
   !      !init mpi
call MPI_init(error)
!      ! get number of processes
call MPI_Comm_size(MPI_COMM_WORLD,numproc,error)
!      ! get individual process id
call MPI_Comm_rank(MPI_COMM_WORLD,id,error)
   
   call directory(id)
   call reader1

   iseed = -8743432
   call readfile('/home/lewis/Desktop/Monte-Carlo-Radiation/data/target.txt', target_a)

   allocate(source(3),test(3))
   
   !set random data for source
   do i = 1,size(source)
      source(i) = (1.d-3 - 1.d-5)*ran2(iseed)+(1.d-5)
   end do
   
!write source out
open(16, file='source.dat')
write(16,*) source
close(16)

   fitval = fitness(source, target_a, error,numproc,id)
   i=0
   
   do while(.True.)
      
      i = i + 1
      
      if(allocated(m))then
         m = mutate(source)
      else
         allocate(m(size(source)))
         m = mutate(source)
      end if
      
      fitval_m = fitness(m, target_a, error,numproc,id)
      
      if(fitval_m .lt. fitval)then
         fitval = fitval_m
         source = m
         print*, i, fitval_m
      end if
      
      if(fitval .lt. 10.d0)then
         open(12,file='values.dat')
         do j=1,size(source)
            write(12,*) source(j)
         end do
         exit
      end if
   
   end do
   call MPI_Finalize(error)
   Contains 
      
      double precision function fitness(concs, tar, error,numproc,id)
         use main
         implicit none
         double precision, dimension(:), intent(IN)  :: tar, concs
         double precision, allocatable, dimension(:) :: source
         integer                                     :: i, error,numproc,id
         !do monty
         call mcpolar(concs, error,numproc,id)
         !read in source
call readfile('/home/lewis/Desktop/Monte-Carlo-Radiation/data/fluro_out.dat', source)
         
         fitness = 0
         do i = 1, size(tar)
            fitness = fitness + (tar(i) - source(i))**2.
         end do
         print*,fitness
     end function fitness
      
     function mutate(source)
     
      use random
      
      double precision, dimension(:), intent(in) :: source
      double precision                           ::  mutate(size(source))
      real                                       :: ran2
      integer                                    :: pos
      
      pos =  int(size(source) * ran2(iseed)) + 1
      mutate = source
      mutate(pos) = (1.d-3 - 1.d-5)*ran2(iseed)+(1.d-5)

      end function
           
 subroutine readfile(filename, array)

   double precision, allocatable, dimension(:), intent(inout) :: array
   double precision                                           :: a
   character(*),                                intent(in)    :: filename
   integer                                                    :: cnt, io, i
   
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
