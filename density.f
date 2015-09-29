      subroutine density(x,y,z,rho)

      implicit none

      real x,y,z,rho

      real w,w2,r,r2

c***** calculate some distances for use in setting up density 
c***** structure. Note that distances are in units of xmax, ymax, and zmax 
c***** as called from the loop over cells in gridset.f
      w2=x*x+y*y
      w=sqrt(w2)
      r2=w2+z*z
      r=sqrt(r2)

c***** Set up uniform density cube within the grid
      if(r.gt.1.) then
        rho=1.
      else
        rho=1.
      endif


      return
      end

