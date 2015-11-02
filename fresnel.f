      subroutine fresnel(n1,n2,cost,sint,sflag,tflag,iseed,
     +      reflc,xcell,ycell,cnt,trans,weight)
      
      implicit none
      real n1,n2,cost,tir,cost2,sint,f1,f2,ran2,weight
      integer iseed,cnt,xcell,ycell
      logical sflag,tflag
      real reflc(1:cnt,1:cnt),trans(1:cnt,1:cnt)
     
     
      if(sflag.eqv..TRUE.)then
            if(n1.eq.n2)then!equal refrative indices
                  tir=0.
            elseif(abs(cost).ge.1.-1.E-8)then!cost straight down
                  tir=(n1-n2)**2/(n2+n1)**2
            elseif(abs(cost).lt.1.E-6)then!oblique angle
                  tir=1.0
            else
!                  cost=abs(cost)
                  sint=sqrt(1.-cost**2)
                  if(n2*sint/n1.gt.1.)print *,'shit'
                  cost2=sqrt(1.-(n1*sint/n2)**2)

                  f1=abs((n1*cost-n2*cost2)/(n1*cost+n2*cost2))**2
                  f2=abs((n1*cost2-n2*cost)/(n1*cost2+n2*cost))**2
                  tir=0.5*(f1+f2)
            end if
!            print *, tir
            if(ran2(iseed).lt.tir)then
                  tflag=.TRUE.
                       reflc(xcell,ycell)=reflc(xcell,ycell)+weight
            else
                  trans(xcell,ycell)=trans(xcell,ycell)+weight
                  weight=0.

      end if
            end if
      return       
      end
            
     
