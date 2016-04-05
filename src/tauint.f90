MODULE tauint


   implicit none
   save
   
CONTAINS

   recursive subroutine tauint2(n1,n2,xcell,ycell,zcell,&
                     tflag,iseed,delta,sflag,weight,ddx,ddy)

   use constants,   only : PI,nxg,nyg,nzg,xmax,ymax,zmax,OFFSET
   use photon_vars, only : xp,yp,zp,nxp,nyp,nzp
   use iarray,      only : jmean,xface,yface,zface,rhokap!,trans,fluroexit
   use opt_prop,    only : kappa
   use noisey_mod

   implicit none

   integer iseed,xcell,ycell,zcell
   DOUBLE PRECISION weight,ddx,ddy
   DOUBLE PRECISION n1,n2
   logical tflag,sflag
   real ran2

   integer celli,cellj,cellk
   DOUBLE PRECISION tau,taurun,taucell,d,d1,dcell,xcur,ycur,zcur,dsx,dsy,dsz
   DOUBLE PRECISION dx,dy,dz,smax,delta

!**** generate random optical depth tau
tau=-alog(ran2(iseed))

!***** set the cumulative distance and optical depth (d and taurun) 
!***** along the photon path to zero.  set the current photon coordinates.
!***** note that the origin of the (xcur,ycur,zcur) system is at the 
!***** bottom corner of the grid.
   taurun=0.
   d=0.

   xcur=xp+xmax
   ycur=yp+ymax
   zcur=zp+zmax

   celli=xcell
   cellj=ycell
   cellk=zcell

!***** calculate smax -- maximum distance photon can travel
   if(nxp.gt.0.) then
      dsx=(2.*xmax-xcur)/nxp
   elseif(nxp.lt.0.) then
      dsx=-xcur/nxp
   elseif(nxp.eq.0.) then
      dsx=1.e2*xmax
   endif

   if(nyp.gt.0.) then
      dsy=(2.*ymax-ycur)/nyp
   elseif(nyp.lt.0.) then
      dsy=-ycur/nyp
   elseif(nyp.eq.0.) then
      dsy=1.e2*ymax
   endif

   if(nzp.gt.0.) then
      dsz=(2.*zmax-zcur)/nzp
   elseif(nzp.lt.0.) then
      dsz=-zcur/nzp
   elseif(nzp.eq.0.) then
      dsz=1.e2*zmax
   endif

   smax=min(dsx,dsy,dsz)
   if(smax.lt.delta) then
      tflag=.TRUE.
      return
   endif
 
!***** integrate through grid
   do while((taurun.lt.tau).and.(d.lt.(.999*smax)))

!***** find distance to next x, y, and z cell walls.  
!***** note that dx is not the x-distance, but the actual distance along 
!*****the direction of travel to the next x-face, and likewise for dy and dz.
   if(nxp.gt.0.) then
   dx=(xface(celli+1)-xcur)/nxp
   if(dx.lt.delta) then
      xcur=xface(celli+1)
      celli=celli+1
      dx=(xface(celli+1)-xcur)/nxp
   endif
      elseif(nxp.lt.0.) then
   dx=(xface(celli)-xcur)/nxp
   if(dx.lt.delta) then
      xcur=xface(celli)
      dx=(xface(celli-1)-xcur)/nxp
      celli=celli-1
   endif
      elseif(nxp.eq.0.) then
   dx=1.e2*xmax
      endif

      if(nyp.gt.0.) then
   dy=(yface(cellj+1)-ycur)/nyp
   if(dy.lt.delta) then
      ycur=yface(cellj+1)
      cellj=cellj+1
      dy=(yface(cellj+1)-ycur)/nyp
   endif
      elseif(nyp.lt.0.) then
   dy=(yface(cellj)-ycur)/nyp
   if(dy.lt.delta) then
      ycur=yface(cellj)
      dy=(yface(cellj-1)-ycur)/nyp
      cellj=cellj-1
   endif
      elseif(nyp.eq.0.) then
   dy=1.e2*ymax
      endif

      if(nzp.gt.0.) then
   dz=(zface(cellk+1)-zcur)/nzp
   if(dz.lt.delta) then
      zcur=zface(cellk+1)
      cellk=cellk+1
      dz=(zface(cellk+1)-zcur)/nzp
   endif
      elseif(nzp.lt.0.) then
   dz=(zface(cellk)-zcur)/nzp
   if(dz.lt.delta) then
      zcur=zface(cellk)
      dz=(zface(cellk-1)-zcur)/nzp
      cellk=cellk-1
   endif
      elseif(nzp.eq.0.) then
   dz=1.e2*zmax
      endif

!   if((celli.gt.nxg).or.(celli.lt.0)) print*,'celli',celli,xcur
!   if((cellj.gt.nyg).or.(cellj.lt.0)) print*,'cellj',cellj,ycur
!   if((cellk.gt.nyg).or.(cellk.lt.0)) print*,'cellk',cellk,zcur


!***** distances are only zero if photon is on cell wall.  if it is 
!***** on cell wall then set to arbitrary large distance, since we will
!***** in fact hit another wall

!******* THIS USED TO BE (dx.eq.0.)

   if( (dx.eq.0.) .or. ((abs(dx)).lt.(delta)) ) dx=1.e2*xmax
   if( (dy.eq.0.) .or. ((abs(dy)).lt.(delta)) ) dy=1.e2*ymax
   if( (dz.eq.0.) .or. ((abs(dz)).lt.(delta)) ) dz=1.e2*zmax

!***** find distance to next cell wall -- minimum of dx, dy, and dz
   dcell=min(dx,dy,dz)
   if(dcell.le.0.) then
   print *,'tauint2: dcell < 0',dx,dy,dz,nxp,nyp,nzp
   print *,xcur,ycur,zcur,celli,cellj,cellk
      endif
   if(dx.lt.0.) dcell=min(dy,dz)
   if(dy.lt.0.) dcell=min(dx,dz)
   if(dz.lt.0.) dcell=min(dx,dy)

!***** optical depth to next cell wall is 
!***** taucell= (distance to cell)*(opacity of current cell)
   taucell=dcell*rhokap(celli,cellj,cellk,1)

!***** if taurun+taucell>tau then scatter at distance d+d1.  
!***** update photon position and cell.  
!***** if taurun+taucell<tau then photon moves distance dcell 
!***** (i.e. ends up on next cell wall) and update photon position
!***** and cell.
   if((taurun+taucell).ge.tau) then
   d1=(tau-taurun)rhokap(celli,cellj,cellk,1)
   d=d+d1
   jmean(celli,cellj,cellk,1)=jmean(celli,cellj,cellk,1)+d1
   taurun=taurun+taucell
   xcur=xcur+d1*nxp
   ycur=ycur+d1*nyp
   zcur=zcur+d1*nzp

!*************** Linear Grid ************************
   celli=int(nxg*xcur/(2.*xmax))+1
   cellj=int(nyg*ycur/(2.*ymax))+1
   cellk=int(nzg*zcur/(2.*zmax))+1
!****************************************************

   else

   d=d+dcell
   jmean(celli,cellj,cellk,1)=jmean(celli,cellj,cellk,1)+dcell
   taurun=taurun+taucell
   xcur=xcur+dcell*nxp
   ycur=ycur+dcell*nyp
   zcur=zcur+dcell*nzp

!*************** Linear Grid ************************
   celli=int(nxg*xcur/(2.*xmax))+1
   cellj=int(nyg*ycur/(2.*ymax))+1
   cellk=int(nzg*zcur/(2.*zmax))+1
!****************************************************
    endif

   end do

!***** calculate photon final position.  if it escapes envelope then
!***** set tflag=.TRUE.  if photon doesn't escape leave tflag=.FALSE. and update 
!***** photon position.
   if((d.ge.(.999*smax))) then

      
      if(zcur.gt.(2.*zmax)-OFFSET)then
         call fresnel(n1,n2,sflag,tflag,iseed,ddx,ddy,weight,xcur,ycur)
      else
            tflag=.TRUE.
      end if

   else

      xp=xp+d*nxp
      yp=yp+d*nyp
      zp=zp+d*nzp
      xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
      ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
      zcell=int(nzg*(zp+zmax)/(2.*zmax))+1

   endif

   end subroutine tauint2
   
   
!*************************************************************************   


   subroutine fresnel(n1,n2,sflag,tflag,iseed,ddx,ddy,weight,xcur,ycur)
        
   use constants, only : pi,cbinsnum
   use photon_vars, only : nxp,nyp,nzp,cost,sint,cosp,sinp
!   use iarray, only : trans

   implicit none
        

   DOUBLE PRECISION n1,n2,n,tir,cost2,f1,f2,xcur,ycur
   DOUBLE PRECISION costt,crit,ran,weight,ddx,ddy
   real ran2
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
        !call tauint2(n1,n2,xcell,ycell,zcell,&
         !      tflag,iseed,delta,sflag,weight,ddx,ddy)
   else
   !photon transmitted
   if(nzp.gt.0.)then
   ix=floor(xcur/ddx)
   iy=floor(ycur/ddy)
   if(ix.gt.cbinsnum)ix=cbinsnum
   if(ix.lt.1)ix=1
   if(iy.gt.cbinsnum)iy=cbinsnum
   if(iy.lt.1)iy=1            
!   trans(ix,iy)=trans(ix,iy)+weight
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

   end subroutine fresnel   
end MODULE tauint
