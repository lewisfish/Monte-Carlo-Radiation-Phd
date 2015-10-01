      subroutine density(x,y,z,rho,kflag,kappa1,kappa2)

      implicit none

      real x,y,z,rho,kappa1,kappa2
      integer kflag
      real w,w2,r,r2

c***** calculate some distances for use in setting up density 
c***** structure. Note that distances are in units of xmax, ymax, and zmax 
c***** as called from the loop over cells in gridset.f
!      w2=x*x+y*y
!      w=sqrt(w2)
      r2=x*x+y*y+z*z
      r=sqrt(r2)

c***** Set up optically diffrent sphere within the grid
      if(kflag.eq.1)then
            if(r.gt..4) then
                  rho=kappa1
            else
                  rho=kappa2
            endif
      else
            if(r.gt..4) then
                  rho=kappa2
            else
                  rho=kappa2
            endif
      end if
      return
      end

