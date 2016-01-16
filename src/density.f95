MODULE density_mod

implicit none
save

CONTAINS
   subroutine density(x,y,z,rho,kappa,cur,xmax,ymax,zmax)

   implicit none

   integer cur
   real x,y,z,rho,kappa(cur)
   real xmax,ymax,zmax

!***** calculate some distances for use in setting up density 
!***** structure. Note that distances are in units of xmax, ymax, and zmax 
!***** as called from the loop over cells in gridset.f



!***** Set up optically diffrent sphere within the grid
   if((x.gt.-2.*xmax .and. x.lt.2.*xmax).and.(y.gt.-2.*ymax .and. y.lt.2.*ymax).and.(z.gt.-2.*zmax .and. z.lt.2.*zmax))then
         rho=kappa(1)
   else  
         rho=0.
   end if

   return
   end
end MODULE density_mod

