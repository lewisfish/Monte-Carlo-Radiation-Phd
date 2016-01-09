subroutine facecheck(face,celli,cellj,cellk,iseed,
     +                    nxp,nyp,nzp)

implicit none

integer,intent(IN) :: face(6),celli,cellj,cellk
integer :: temp(6),loctemp,iseed
integer:: location
double precision :: angle
real ran2
double precision,intent(INOUT) :: nxp,nyp,nzp

!check to see if photon clode to crystal
if((celli.gt.face(1)-5.and.celli.lt.face(2)+5).and.
     +(cellj.gt.face(3)-5.and.cellj.lt.face(4)+5).and.
     +(cellk.gt.face(5).and.cellk.lt.face(6)-5))then
     
!calculate the distance between photon cell and crystal cell
! so we can work out which face is closest

temp(1)=abs(face(1)-celli)
temp(2)=abs(face(2)-celli)
temp(3)=abs(face(3)-cellj)      
temp(4)=abs(face(4)-cellj)
temp(5)=abs(face(5)-cellk)
temp(6)=abs(face(6)-cellk)

!check to see if photon on crystal edge
if(count(temp.eq.0).gt.1)then
     location=minloc(temp,DIM=1)
     if(ran2(iseed).gt..5)then
     temp(location)=10
     end if
end if

if(minval(temp,DIM=1).eq.0)then
      location=minloc(temp,DIM=1)
      call facenormal(nxp,nyp,nzp,location,angle)
end if

end if
end subroutine facecheck
