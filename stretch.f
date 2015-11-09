      subroutine stretch(sfact,nxp,nyp,nzp,p,stretchflag)
      
      implicit none
      
      double precision,intent(OUT) :: sfact
      double precision :: mu,mag
      double precision,intent(IN) :: nxp,nyp,nzp,p
      logical,intent(out) :: stretchflag
      
      !calculate the cosine of the angle between the prefered 
      !direction and the direction the photon is going in
      mag=sqrt(nxp**2+nyp**2+nzp**2)
      mu=cos((-1.*nzp)/(mag))
      
      !calculate the stretch factor
      sfact=p*mu
      stretchflag=.true.
      end subroutine stretch
