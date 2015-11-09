      subroutine density(x,y,z,rho,kflag,kappa,cur,xlow,xhi,ylow,
     +                  yhi,zlow,zhi)

      implicit none

      integer kflag,cur
      double precision x,y,z,rho,kappa(cur)
      double precision xlow,xhi,ylow,yhi,zlow,zhi
      
c***** calculate some distances for use in setting up density 
c***** structure. Note that distances are in units of xmax, ymax, and zmax 
c***** as called from the loop over cells in gridset.f



c***** Set up optically diffrent sphere within the grid
            if((y.lt.yhi.and.y.gt.ylow).and.(z.lt.zhi.and.z.gt.zlow)
     +      .and.(x.lt.xhi.and.x.gt.xlow))then
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

