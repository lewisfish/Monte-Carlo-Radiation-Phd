program error
   implicit none

   real,allocatable :: array(:,:)
   real :: a, b
   integer :: io, i, nzg
   
   open(78,file='fit_params.txt', status = 'old' , IOSTAT = io )
      read(78,*) a, b
   close(78)
   call readfile_array2D("testff.dat", array, 0, 2)
   nzg=size(array) / 2
   open(45, file = "error.dat")
   do i=1, nzg
!      write(45,*) array(i,1) ,abs(1. - (array(i,2)/(.771*exp(-.994*array(i,1)/.29))))
      write(45,*) array(i,1), abs(1.-(a*exp(-b*array(i,1)/.29))/(.771*exp(-.994*array(i,1)/.29)))
   end do
   
   contains
   
      subroutine readfile_array2D(filename, array, flag, colsize)
      !
      ! Reads a file to get its length and allocates an array to store the data and reads it in.
      !
      ! subroutine takes filename, the array for dat to be read into, a flag to toggle
      ! between using square array and other shapes, colsize is an optional argument that
      ! specifies the size of the 2nd dimension of the array
      !
      real, allocatable, dimension(:,:), intent(inout) :: array
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
end program
