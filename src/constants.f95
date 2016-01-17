MODULE constants
!
! Module containing constants:
!         PI,TWOPI, the number of grid elements in each direction n%g,
!         the various bin # params,
!         and the vars for the filepaths.
!
implicit none
save

integer, parameter :: nxg=200,nyg=200,nzg=200,Nbins=401,cbinsnum=200
real, parameter    :: PI = 3.141592,TWOPI=6.283185
character(len=255) :: cwd,homedir,fileplace,resdir

end MODULE constants
