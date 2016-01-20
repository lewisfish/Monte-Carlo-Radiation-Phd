MODULE binning_mod

implicit none
save

CONTAINS
   subroutine binning(ddr,zcur,ddz,absorb)

   use constants, only   : cbinsnum
   use photon_vars, only : xp,yp
   use iarray, only      : deposit
   
   implicit none

   integer :: binz,binr
   real    :: zcur,ddr,ddz,absorb

   binr = floor(sqrt(xp**2 + yp**2)/ddr)
   binz = floor(zcur/ddz)

   if(binz.lt.1..or.binr.lt.1.)then
   !do nothing
   elseif(binz.gt.cbinsnum.or.binr.gt.cbinsnum)then
   !do nothing
   else
      deposit(binr,binz) = deposit(binr,binz) + absorb
   end if

   end subroutine binning
end MODULE binning_mod
