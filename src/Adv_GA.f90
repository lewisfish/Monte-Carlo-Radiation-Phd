MODULE random

   integer                       :: iseed, elite, gensize
   double precision, allocatable :: target_a(:)

end module

program simple_GA

   use mpi
   use random
   use subs
   use reader_mod
   use main
   implicit none

   double precision              :: cf, minfit
   real                          :: ran2,start,finish
   double precision, allocatable :: genepool(:, :)
   double precision, allocatable :: fit_array(:), child(:), parent1(:), parent2(:)
   integer                       :: i, j, maxl, st, stall_counter, termination_val
   character(len=255)            :: filen 
   integer                       :: id,error,numproc
   
! init mpi
call MPI_init(error)
! get number of processes
call MPI_Comm_size(MPI_COMM_WORLD,numproc,error)
! get individual process id
call MPI_Comm_rank(MPI_COMM_WORLD,id,error)
   
   call cpu_time(start)
   termination_val = 500
   if(id.eq.0)print*,'Starting GA on',numproc,'cores'
   filen = '/home/lewis/Desktop/Ga-simple-monte/old/data/'
   
   open(95,file=trim(filen)//'fitvals.dat',status='replace')
   write(95,*)''
   close(95)

   call directory(id)
   call reader1

   
   !seed the random number generator
   iseed = -8743432
   
   !read in target data and allocate parent and child sizes
   call readfile(trim(filen)//'target.txt', target_a)
   st = 3!size(target_a)
   allocate(child(st), parent1(st), parent2(st))

   !# of parents in genepool
   gensize = 10
   elite = 3+int(0.01*gensize)
   allocate(genepool(st+1, gensize))
   
   !init genepool woth random data
   if(id.eq.0)then
   do i = 1, gensize
      do j = 1, st
         genepool(j, i) = (1.d-3 - 1.d-5)*ran2(iseed)+(1.d-5)
      end do
      genepool(st+1,i) = fitness(genepool(:st+1, i) , target_a, .false.)
      if(id.eq.0)print*,'added ',i,genepool(st+1,i)
   end do
   call MPI_BCAST(genepool,(st+1)*gensize,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error)
   print*,'Populated Genepool...'
   end if
   call MPI_Barrier(MPI_COMM_WORLD,error)
   i=0
   minfit=1.d30
   stall_counter = 0

   do while(.True.)
      
      i = i + 1   !generation number
   call MPI_Barrier(MPI_COMM_WORLD,error)
      !if the smallest fitness value in the genepool is less than threshold
      !then exit
      if(minval(genepool(st+1,:)) .lt. 0.01 .or. &
         stall_counter .ge. termination_val)then
         !target reached
         call MPI_Barrier(MPI_COMM_WORLD,error)
         if(id.eq.0)then
            open(12,file='output.dat')
            print*,i,minval(genepool(st+1,:)),stall_counter
            do j=1, st
               write(12,*) genepool(j,int(minloc(genepool(st+1,:), 1)))
            end do
         end if
         call MPI_Barrier(MPI_COMM_WORLD,error)
         exit
      end if
           
      call sort(genepool)
      
      call selection(genepool)
 
      if(minval(genepool(st+1,:)) .lt. minfit)then
         minfit = minval(genepool(st+1,:))
         stall_counter = 0
         call MPI_Barrier(MPI_COMM_WORLD,error)
         if(id.eq.0)then
            print*,genepool(:, int(minloc(genepool(st+1,:), 1)))
            open(12,file='output.dat', status='old')
            do j=1, st
               write(12,*) genepool(j,int(minloc(genepool(st+1,:), 1)))
            end do
            close(12)
         end if
      else
         stall_counter = stall_counter + 1
      end if
      
      if(id.eq.0)then
         print*,i,minfit

         open(95,file='fitvals.dat',access='append')
         write(95,*) i,minfit
         close(95)    
      end if
      call MPI_Barrier(MPI_COMM_WORLD,error)
   end do
   if(id.eq.0)then
      call cpu_time(finish)
      print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
   end if
   call MPI_Finalize(error)
   Contains 
      
      integer function random_parent()
      !
      !  get a random parent form genepool
      !
         use random
      
         implicit none
         
         integer :: rand
         
         rand = int(ran2(iseed) * ran2(iseed) * real(gensize - 1))+1
         random_parent = rand
         
      end function random_parent
      
      double precision function fitness(concs, tar, pflag)
      !
      !  calculate fitness
      !
         double precision, dimension(:), intent(IN)  :: tar, concs
         double precision, allocatable, dimension(:) :: source
         integer :: i
         logical :: pflag
         
         call MPI_BCAST(concs,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,error)
         call MPI_Barrier(MPI_COMM_WORLD,error)
         call mcpolar(concs,error,numproc,id, pflag)

         call readfile(trim(filen)//'fluro_out.dat', source)
         
         fitness = 0
         do i = 1, size(tar)
            fitness = fitness + (tar(i) - source(i))**2.
         end do
         
     end function fitness
      
     function mutate(parent1, parent2)
     !
     !   create child and mutate it
     !
      use random, only : iseed, target_a
      
      double precision, dimension(:), intent(in) :: parent1, parent2
      double precision                           :: mutate(size(parent2)+1)
      integer                        :: start, stopp, tmp, pos
      real                           :: ran2
      print*,size(parent1),size(parent2),size(mutate)
      mutate = parent1
      !mix dna
      start = int(size(parent1) * ran2(iseed)) + 1
      stopp = int(size(parent2) * ran2(iseed)) + 1
      if(start .gt. stopp)then
         tmp = stopp
         stopp = start
         start = tmp
      end if
      mutate(start:stopp) = parent2(start:stopp)
      
      !mutate child dna
      pos = int(size(mutate) * ran2(iseed)) + 1
      mutate(pos) = (1.d-3 - 1.d-5)*ran2(iseed)+(1.d-5)

      !calc fitness of new child
      mutate(size(mutate)) = fitness(mutate(:size(mutate)-1), target_a,.false.)
      end function
           
     subroutine selection(a)          
   
   use random, only : iseed, gensize, elite
   
   implicit none
           
   double precision, dimension(:,:) :: a
   double precision                 :: tmp(size(a, 1), size(a, 2))
   integer                          :: parent1, parent2, i
   
   tmp=0.
   !elite selection
   tmp(:,:elite) = a(:,:elite)

   do i = elite+1, gensize
      !pick parents randomly
      parent1 = int(ran2(iseed)*(size(a,2)-elite)) + 1
      parent2 = int(ran2(iseed)*(size(a,2)-elite)) + 1
      !add child to new generation
      tmp(:,i) = mutate(a(:size(a,1)-1,parent1), a(:size(a,1)-1,parent2))
   end do
   !sort new generation
   call sort(tmp)
   genepool = tmp     

 end subroutine selection   
 
        
 subroutine sort(a)
 !
 ! sort 2d array by last column
 !
   implicit none
   
    double precision, dimension(:,:), Intent(INOUT) :: a
    integer                                         :: i, minIndex
    double precision                                :: temp(size(a,1))
 
    do i = 1, size(a, 2)-1
       minIndex = MINLOC(a(st+1, i:), 1) + i - 1
       if(a(st+1, i) .gt. a(st+1, minIndex))then
          temp = a(:, i)
          a(:, i) = a(:, minIndex)
          a(:, minIndex) = temp
       end if
    end do

 end subroutine sort      
           
           
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
