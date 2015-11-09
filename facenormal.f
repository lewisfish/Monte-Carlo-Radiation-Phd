      subroutine facenormal(nxp,nyp,nzp,location,angle)
      
      implicit none
      
      double precision,intent(INOUT) :: nxp,nyp,nzp
      integer,intent(IN) :: location
      double precision :: angle
      
      !calculate the angle of inc(N.V)
      if(location.eq.1)then
            if(nxp.gt.0.)then
                  angle=nxp
            else
                  angle=-nxp
            end if
      elseif(location.eq.2)then
            if(nyp.gt.0.)then
                  angle=-nxp
            else
                  angle=nxp
            end if
      elseif(location.eq.3)then
            if(nyp.gt.0.)then
                  angle=nyp
            else
                  angle=-nyp
            end if
      elseif(location.eq.4)then
            if(nyp.gt.0.)then
                  angle=-nyp
            else
                  angle=nyp
            end if
      elseif(location.eq.5)then
            if(nzp.gt.0.)then
                  angle=nzp
            else
                  angle=-nzp
            end if
      elseif(location.eq.6)then
            if(nzp.gt.0.)then
                  angle=-nzp
            else
                  angle=nzp
            end if
      end if
      
       if(n1.eq.n2)then!equal refrative indices
            tir=0.
      elseif(abs(angle).ge.1.*.999)then!cost straight into
            tir=(n1-n2)**2/(n2+n1)**2
      elseif(abs(angle).lt.1.E-6)then!oblique angle
            tir=1.0
      else
            sint=sqrt(1.-angle**2)
            if(n1*sint/n2.gt.1.)print *,'shit',sint,n1,n2
            cost2=sqrt(1.-(n1*sint/n2)**2)

            f1=abs((n1*angle-n2*cost2)/(n1*angle+n2*cost2))**2
            f2=abs((n1*cost2-n2*angle)/(n1*cost2+n2*angle))**2
            tir=0.5*(f1+f2)
      end if
      
      if(ran2(iseed).lt.tir)then
            !light reflected
                  angle=-angle
            else
                  !light transmitted
                  !change opt properties
            end if
      
      
      
      end subroutine facenormal
