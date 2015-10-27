      subroutine density(x,y,z,rho,kflag,kappa,cur)

      implicit none

      integer kflag,i,j,cur
      real x,y,z,rho,kappa(cur)
      real r,r2
      
c***** calculate some distances for use in setting up density 
c***** structure. Note that distances are in units of xmax, ymax, and zmax 
c***** as called from the loop over cells in gridset.f


!      r2=(x+pos(i,j))**2+y*y+(z+pos(i,j))**2
!      r=sqrt(r2)

c***** Set up optically diffrent sphere within the grid
            if((y.lt..3.and.y.gt.-.3).and.(z.lt.-0.1.and.z.gt.-0.4)
     +      .and.(x.lt..5.and.x.gt.-.5))then
                  if(kflag.eq.1)rho=kappa(2)*4.56
                  if(kflag.eq.2)rho=kappa(4)*4.56
                  if(kflag.eq.3)rho=kappa(6)*4.56
                  if(kflag.eq.4)rho=kappa(8)*4.56
            else
            
                  if(kflag.eq.1)rho=kappa(1)*1.113
                  if(kflag.eq.2)rho=kappa(3)*1.113
                  if(kflag.eq.3)rho=kappa(5)*1.113
                  if(kflag.eq.4)rho=kappa(7)*1.113
            end if

      return
      end

