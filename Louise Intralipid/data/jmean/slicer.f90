program slicer


implicit none

integer :: i,j,io,ndim
integer,parameter :: nxg=100, nyg=100, nzg=100
DOUBLE PRECISION  :: outer(nxg,nyg,nzg,1)
DOUBLE PRECISION :: jmean(nxg+3,nyg+3,nzg+3,4)
character(len=255) :: s, filepath


!call get_command_argument(1, filepath)
!call get_command_argument(2, s)
!read(s,*) ndim

!!!**** Read in jmean grid as unformatted array
inquire(iolength=i)jmean
print*,i
open(10,file='jmean.dat',form='unformatted',access='direct',recl=((ndim+3)*(ndim+3)*(ndim+3))*8*4,iostat=io)
   if(io.eq.0)then
     read(10,rec=1) jmean
   end if
close(10)

     open(34,file='blue.dat')
     do i=1,ndim
        write(34,*) (jmean(i,50,j,1),j=1,ndim)
     end do
     

open(10,file=trim(filepath),form='unformatted',access='direct',recl=((ndim)*(ndim)*(ndim))*8)      
      write(10,rec=1) jmean(1:nxg,1:nyg,1:nzg,1:1)
close(10)

open(23,file=trim(filepath),form='unformatted',access='direct',recl=((ndim)*(ndim)*(ndim))*8) 
   read(23,rec=1) outer
   inquire(iolength=i) outer
     open(34,file='red.dat')
     do i=1,ndim
        write(34,*) (outer(i,50,j,1),j=1,ndim)
     end do
end
