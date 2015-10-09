      subroutine density(x,y,z,rho,kflag,kappa1,kappa2)

      implicit none

      real x,y,z,rho,kappa1,kappa2
      integer kflag,i,j
      real r,r2,pos(1:2,1:2)

c***** calculate some distances for use in setting up density 
c***** structure. Note that distances are in units of xmax, ymax, and zmax 
c***** as called from the loop over cells in gridset.f

      pos(1,1)=-.5
      pos(1,2)=-.5
      pos(2,1)=0.
      pos(2,2)=.6
!      do i=1,2
!       do j=1,2
!      r2=(x+pos(i,j))**2+y*y+(z+pos(i,j))**2
!      r=sqrt(r2)

c***** Set up optically diffrent sphere within the grid
      if(kflag.eq.1)then
            if((y.lt..2.and.y.gt.-.2).and.(z.lt..2.and.z.gt.-.2))then
                  rho=kappa2
            else
                  rho=kappa1
            end if
      else
            if((y.lt..2.and.y.gt.-.2).and.(z.lt..2.and.z.gt.-.2))then
                  rho=kappa2
            else
                  rho=kappa2
            endif
      end if
!      end do
!      end do
      return
      end

