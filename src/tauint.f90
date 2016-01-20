MODULE tauint


   implicit none
   save
   
CONTAINS

   subroutine tauint2(xmax,ymax,zmax,n1,n2,xcell &
            ,ycell,zcell,tflag,iseed,delta,sflag,weight,ddx,ddy,cnt)

   use constants,   only : PI,nxg,nyg,nzg,cbinsnum
   use photon_vars, only : xp,yp,zp,nxp,nyp,nzp,cost,sint,cosp,sinp
   use iarray,      only : noise,jmean,xface,yface,zface,rhokap,trans,fluroexit
   use opt_prop,    only : wave
   use fresnel_mod
   use noisey_mod

   implicit none

   integer iseed,xcell,ycell,zcell,cnt
   real xmax,ymax,zmax,weight,ddx,ddy
   real ran2,n1,n2
   logical tflag,sflag

   integer celli,cellj,cellk
   real tau,taurun,taucell,d,d1,dcell,xcur,ycur,zcur,dsx,dsy,dsz
   real dx,dy,dz,smax,delta

!***** tflag=.FALSE. means photon is in envelope
   tflag=.FALSE.

!**** generate random optical depth tau
tau=-alog(ran2(iseed))

!***** set the cumulative distance and optical depth (d and taurun) 
!***** along the photon path to zero.  set the current photon coordinates.
!***** note that the origin of the (xcur,ycur,zcur) system is at the 
!***** bottom corner of the grid.
   taurun=0.
   d=0.
   100   continue
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

   if((celli.gt.nxg).or.(celli.lt.0)) print*,'celli',celli,xcur
   if((cellj.gt.nyg).or.(cellj.lt.0)) print*,'cellj',cellj,ycur
   if((cellk.gt.nyg).or.(cellk.lt.0)) print*,'cellk',cellk,zcur


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
   d1=(tau-taurun)/rhokap(celli,cellj,cellk,1)
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
   
   if(zcur.gt.2.*zmax*.999)then
      fluroexit(int(wave))=fluroexit(int(wave))+1
   end if
!   call noisey(xcell,ycell,cnt)
!   call fresnel(n1,n2,sflag,tflag,iseed,ddx,ddy,weight,xcur,ycur)
!      if(nzp.gt.0.)then
!!                  if(tau-taurun.gt.0.)then
!!                        goto 100
!!                  end if
!            tflag=.TRUE.
!      end if
!   else
      tflag=.TRUE.
!   end if
   else
      xp=xp+d*nxp
      yp=yp+d*nyp
      zp=zp+d*nzp
      xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
      ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
      zcell=int(nzg*(zp+zmax)/(2.*zmax))+1

   endif

   end subroutine tauint2
end MODULE tauint
