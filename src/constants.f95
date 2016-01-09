MODULE constants
!
! Module containing constants:
!         PI,TWOPI and the number of grid elements in each direction n%g
! 
! also contains (for now (8/1/16) will move them to array module): 
!         face postions,grid opacity and fluence grid (both local and global)
!
implicit none
save

integer,parameter :: nxg=200,nyg=200,nzg=200
real,parameter :: PI = 3.141592,TWOPI=6.283185
character(len=255) :: cwd,homedir,fileplace,resdir

real :: xface(nxg+3),yface(nyg+3),zface(nzg+3)
real :: rhokap(nxg+3,nyg+3,nzg+3,4)
real :: jmean(nxg+3,nyg+3,nzg+3,4)
real :: jmeanGLOBAL(nxg+3,nyg+3,nzg+3,4)

end MODULE constants
