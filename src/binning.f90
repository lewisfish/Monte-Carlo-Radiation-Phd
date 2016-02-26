MODULE binning_mod

implicit none
save

CONTAINS
   subroutine binning(weight)

   use constants, only   : cbinsnum,nxg,nyg,nzg,xmax,ymax,zmax
   use photon_vars, only : xp,yp,zp
   use iarray, only      : deposit
   
   implicit none

   integer :: i, j, k
   real, optional, intent(in) :: weight

   i = int(nxg * (xp+xmax) / (2.*xmax)) + 1
   j = int(nyg * (yp+ymax) / (2.*ymax)) + 1
   k = int(nzg * (zp+zmax) / (2.*zmax)) + 1
   
   if(present(weight))then
      deposit(i, j, k) = deposit(i, j, k) + weight
   else
      deposit(i, j, k) = deposit(i, j, k) + 1
   end if
   end subroutine binning
end MODULE binning_mod
