subroutine fresnel(n1,n2,sflag,tflag,iseed,ddx,ddy,weight,xcur,ycur)
     
use constants, only : pi,cbinsnum
use photon, only : nxp,nyp,nzp,cost,sint,cosp,sinp
use iarray, only : trans

implicit none
     

real n1,n2,n,tir,cost2,f1,f2,xcur,ycur
real ran2,costt,crit,ran,weight,ddx,ddy
integer iseed,ix,iy
logical sflag,tflag
     
     
!swap refractive index if going in to medium
if(sflag.eqv..TRUE.)then
n=0.
n=n1
n1=n2
n2=n

end if
     
crit=asin(n2/n1)
!      print*,crit*180./pi 
if(nzp.lt.0.)then!adjust angle for fresnel calculation when photon heading down
costt=abs(cost)     
else
costt=cost
end if

!      print*,crit*180./pi,costt
if(n1.eq.n2)then!equal refrative indices
tir=0.
elseif(abs(costt).ge.1.*.999)then!cost straight down
tir=(n1-n2)**2/(n2+n1)**2
!      print*,'straight',acos(costt)*180./pi,tir
elseif(acos(costt).gt.crit)then!total internal reflection
tir=1.0
!            print*,'tir',acos(costt)*180./pi,tir
elseif(abs(costt).lt.1.E-6)then!oblique angle
tir=1.0
!            print*,'oblique',acos(costt)*180./pi,tir
else
sint=sqrt(1.-costt**2)
if(n1*sint/n2.gt.1.)print *,'shit1',sint,n1,n2
cost2=sqrt(1.-(n1*sint/n2)**2)

f1=abs((n1*costt-n2*cost2)/(n1*costt+n2*cost2))**2
f2=abs((n1*cost2-n2*costt)/(n1*cost2+n2*costt))**2
tir=0.5*(f1+f2)

!            print*,'other',acos(costt)*180./pi,tir
end if

ran=ran2(iseed)
if(ran.lt.tir)then
!           photon reflected
     cost=-cost
else
!photon transmitted
if(nzp.gt.0.)then
ix=floor(xcur/ddx)
iy=floor(ycur/ddy)
if(ix.gt.cbinsnum)ix=cbinsnum
if(ix.lt.1)ix=1
if(iy.gt.cbinsnum)iy=cbinsnum
if(iy.lt.1)iy=1            
trans(ix,iy)=trans(ix,iy)+weight
      if(sflag.eqv..FALSE.)then
             tflag=.TRUE.
      end if
end if

!            weight=0.
end if

sint=sqrt(1.-cost**2)

nxp=sint*cosp  
nyp=sint*sinp
nzp=cost

!swap index back
if(sflag.eqv..TRUE.)then
n2=n1
n1=n
end if

return       
end

     
