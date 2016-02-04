MODULE binning_mod

implicit none
save

CONTAINS
   subroutine binning(xmax,ymax,zmax,ddz,xcur,ycur,zcur)

   use constants, only   : cbinsnum,nxg,nyg,nzg
   use photon_vars, only : xp,yp,zp
   use iarray, only      : deposit
   
   implicit none

   integer :: celli,cellj,cellk
   real    :: zmax,ddr,ddz,xcur,ycur,zcur,xmax,ymax

   celli=int(nxg*xcur/(2.*xmax))+1
   cellj=int(nyg*ycur/(2.*ymax))+1
   cellk=int(nzg*zcur/(2.*zmax))+1
   
!   if(binz.lt.1..or.binx.lt.1..or.biny.lt.1.)then
!   !do nothing
!   elseif(binz.gt.cbinsnum.or.binx.gt.cbinsnum.or.biny.gt.cbinsnum)then
!   !do nothing
!   else
      deposit(celli,cellj,cellk) = deposit(celli,cellj,cellk) + 1
!   end if

   end subroutine binning
end MODULE binning_mod
