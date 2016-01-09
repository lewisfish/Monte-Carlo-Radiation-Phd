subroutine binning(deposit,ddr,cbinsnum,zcur,ddz,absorb)

use photon,only : xp,yp

implicit none

integer binz,binr,cbinsnum
real  deposit(cbinsnum,cbinsnum)
real zcur,ddr,ddz,absorb

      binr=floor(sqrt(xp**2+yp**2)/ddr)
      binz=floor(zcur/ddz)

if(binz.lt.1..or.binr.lt.1.)then
!do nothing
elseif(binz.gt.cbinsnum.or.binr.gt.cbinsnum)then
!do nothing
else
deposit(binr,binz)=deposit(binr,binz)+absorb
end if
end subroutine binning
