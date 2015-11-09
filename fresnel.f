      subroutine fresnel(n1,n2,cost,sint,sflag,tflag,iseed,
     +      reflc,xcell,ycell,cnt,trans,weight,pi,nxp,nyp,nzp,
     +      cosp,sinp)
      
      implicit none
      double precision n1,n2,cost,tir,cost2,sint,f1,f2
      real ran2
      double precision nyp,nzp,sinp,cosp,weight,pi,nxp
      integer iseed,cnt,xcell,ycell
      logical sflag,tflag,dir
      double precision reflc(1:cnt,1:cnt),trans(1:cnt,1:cnt)
     
     
      if(sflag.eqv..TRUE.)then
      !light coming into media
      sflag=.FALSE.
      dir=.FALSE.
      else
      !light leaving top face
      dir=.TRUE.
      end if
      
      if(n1.eq.n2)then!equal refrative indices
            tir=0.
      elseif(abs(cost).ge.1.*.999)then!cost straight down
            tir=(n1-n2)**2/(n2+n1)**2
      elseif(abs(cost).lt.1.E-6)then!oblique angle
            tir=1.0
      else
            sint=sqrt(1.-cost**2)
            if(n1*sint/n2.gt.1.)print *,'shit',sint,n1,n2
            cost2=sqrt(1.-(n1*sint/n2)**2)

            f1=abs((n1*cost-n2*cost2)/(n1*cost+n2*cost2))**2
            f2=abs((n1*cost2-n2*cost)/(n1*cost2+n2*cost))**2
            tir=0.5*(f1+f2)
      end if

      if(dir.eqv..FALSE.)then
      !light into media
      print*, 'here, shiit'
            if(ran2(iseed).lt.tir)then
            !light reflected
                  tflag=.TRUE.
                  reflc(xcell,ycell)=reflc(xcell,ycell)+weight
                  weight=0.
            else
                  !light transmitted
                  trans(xcell,ycell)=trans(xcell,ycell)+weight
            end if
      else 
      !light out of media           
            if(ran2(iseed).lt.tir)then
                  ! light reflected
                  cost=-cost
            else
                  !light transmitted
                  tflag=.TRUE.
                  weight=0.
            end if
      end if 
      
      nxp=sint*cosp  
      nyp=sint*sinp
      nzp=cost
           
      return       
      end
            
     
