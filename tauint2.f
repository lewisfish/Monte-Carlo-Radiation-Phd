      subroutine tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,forceflag,
     +    flucount,kappa,xface,yface,zface,rhokap,noise,cnt,d,tau,dens,
     +    p,sfact,xcell,ycell,zcell,tflag,iseed,delta,cur,jmean,
     +    stretchflag,sint,cost,sinp,phi,twopi,tauflag,weight,xexit,
     +    yexit,zexit,pi,sflag,n1,n2,trans,reflc,cosp,face)

      implicit none

      include 'grid.txt'     
     
      integer iseed,xcell,ycell,zcell,cnt,cur,i,face(6),location
      logical tflag,tauflag,forceflag,stretchflag,sflag
      double precision xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,phi,dens(8)
      double precision kappa(8),twopi,sint,cost,sinp,xexit,yexit    
      real ran2
      integer celli,cellj,cellk,flucount
      double precision tau,taurun,taucell,d,d1,dcell,xcur,ycur,zcur,p
      double precision sfact,cosp,n1,n2
      double precision dx,dy,dz,smax,delta,noise,dsx,dsy,dsz,weight,pi
      double precision reflc(1:cnt,1:cnt),trans(1:cnt,1:cnt),zexit

c***** tflag=0 means photon is in envelope
      tflag=.FALSE.
c**** generate random optical depth tau
      if(forceflag.eqv..FALSE.)then
      tau=-log(ran2(iseed))
      end if

c***** set the cumulative distance and optical depth (d and taurun) 
c***** along the photon path to zero.  set the current photon coordinates.
c***** note that the origin of the (xcur,ycur,zcur) system is at the 
c***** bottom corner of the grid.
      taurun=0.
      d=0.
      xcur=xp+xmax
      ycur=yp+ymax
      zcur=zp+zmax

      celli=xcell
      cellj=ycell
      cellk=zcell

c***** calculate smax -- maximum distance photon can travel
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
       
c***** integrate through grid
      do while((taurun.lt.tau).and.(d.lt.(.999*smax)))
      

c***** find distance to next x, y, and z cell walls.  
c***** note that dx is not the x-distance, but the actual distance along 
c*****the direction of travel to the next x-face, and likewise for dy and dz.
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
c***** distances are only zero if photon is on cell wall.  if it is 
c***** on cell wall then set to arbitrary large distance, since we will
c***** in fact hit another wall
         if( (dx.eq.0.) .or. ((abs(dx)).lt.(delta)) ) dx=1.e2*xmax
         if( (dy.eq.0.) .or. ((abs(dy)).lt.(delta)) ) dy=1.e2*ymax
         if( (dz.eq.0.) .or. ((abs(dz)).lt.(delta)) ) dz=1.e2*zmax

c***** find distance to next cell wall -- minimum of dx, dy, and dz
         dcell=min(dx,dy,dz)
         if(dcell.le.0.) then
            print *,'tauint2: dcell < 0'
         endif
         if(dx.lt.0.) dcell=min(dy,dz)
         if(dy.lt.0.) dcell=min(dx,dz)
         if(dz.lt.0.) dcell=min(dx,dy)

c***** optical depth to next cell wall is 
c***** taucell= (distance to cell)*(opacity of current cell)
         if(cur.eq.1)then
         call stretch(sfact,nxp,nyp,nzp,p,stretchflag)
         taucell=dcell*rhokap(celli,cellj,cellk,cur)*(1.-sfact)         
         else
         taucell=dcell*rhokap(celli,cellj,cellk,cur)
         end if
c***** if taurun+taucell>tau then scatter at distance d+d1.  
c***** update photon position and cell.  
c***** if taurun+taucell<tau then photon moves distance dcell 
c***** (i.e. ends up on next cell wall) and update photon position
c***** and cell.

         if((taurun+taucell).ge.tau) then
            if(cur.eq.1)then
            call stretch(sfact,nxp,nyp,nzp,p,stretchflag)
            d1=(tau-taurun)/(rhokap(celli,cellj,cellk,cur)*(1.-sfact))       
            else
            d1=(tau-taurun)/rhokap(celli,cellj,cellk,cur)
            end if            
            d=d+d1
            taurun=taurun+taucell
            xcur=xcur+d1*nxp
            ycur=ycur+d1*nyp
            zcur=zcur+d1*nzp
            !path length estimator
            jmean(celli,cellj,cellk,cur)=jmean(celli,cellj,cellk,cur)+d1
c*************** Linear Grid ************************
            celli=int(nxg*xcur/(2.*xmax))+1
            cellj=int(nyg*ycur/(2.*ymax))+1
            cellk=int(nzg*zcur/(2.*zmax))+1
c****************************************************
!            call facecheck(face,celli,cellj,cellk,iseed)

         else

            d=d+dcell
            taurun=taurun+taucell
            xcur=xcur+dcell*nxp
            ycur=ycur+dcell*nyp
            zcur=zcur+dcell*nzp
            !path length estimator 
        jmean(celli,cellj,cellk,cur)=jmean(celli,cellj,cellk,cur)+dcell
c*************** Linear Grid ************************
            celli=int(nxg*xcur/(2.*xmax))+1
            cellj=int(nyg*ycur/(2.*ymax))+1
            cellk=int(nzg*zcur/(2.*zmax))+1
c****************************************************
          endif
          !changes the opt prop based
!          if(rhokap(xcell,ycell,zcell,cur).ne.kappa(cur))then
!            do i=1,8
!                  if(rhokap(xcell,ycell,zcell,cur).eq.kappa(i))then
!                        if(mod(i,2).eq.0)then
!                              cur=i/2
!                        else
!                              cur=CEILING(real(i)/2.)
!                        end if
!                  end if
!            end do
!          end if 
      end do

c***** calculate photon final position.  if it escapes envelope then
c***** set tflag=TRUE.  if photon doesn't escape leave tflag=FALSE and update 
c***** photon position. 
      
      if((d.ge.(.999*smax))) then
!             tflag=.TRUE.
!           if((xcur.gt.2.*xmax*.999).or.(xcur.lt.0.0001))xexit=xexit+1
!           if((zcur.gt.2.*zmax*.999))yexit=yexit+1
           if((zcur.gt.2.*zmax*.999))then
                  call fresnel(n1,n2,cost,sint,sflag,tflag,iseed,
     +      reflc,xcell,ycell,cnt,trans,weight,pi,nxp,nyp,nzp,
     +      cosp,sinp)
                  if(tflag.eqv..TRUE.)then
                  zexit=zexit+1
                  else
                  yexit=yexit+1
                  end if
           end if
             if(zcur.gt.99.*zmax*.999)then
!             call noisey(xcur,ycur,noise,cnt)
             print*,'not meant to be here'
             end if
             tflag=.TRUE.   
      else
         xp=xp+d*nxp
         yp=yp+d*nyp
         zp=zp+d*nzp
         xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
         ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
         zcell=int(nzg*(zp+zmax)/(2.*zmax))+1

       call flurosub(flucount,rhokap,iseed,xcell,ycell,zcell,cur
     +                   ,kappa,nxp,nyp,sint,cost,sinp,phi,
     +                   twopi,tauflag,weight,dens)

      endif

      return
      end
