      subroutine flurosub(flucount,rhokap,iseed,xcell,ycell,zcell,cur
     +                   ,kappa,nxp,nyp,sint,cost,sinp,phi,
     +                   twopi,tauflag,weight,dens)
      
      implicit none
      
      include 'grid.txt'
    !subroutine to check if photon is in fluro material, then release a fluro photon if appropriate
      integer xcell,ycell,zcell,cur,iseed,flucount
      logical tauflag
      real ran2
      double precision chance,nxp,nyp,nzp,cost,sint,sinp,phi
      double precision twopi,cosp,ran,kappa(8),weight,dens(8)
      
      chance=.50
      ran=ran2(iseed)

      if(rhokap(xcell,ycell,zcell,cur).ne.kappa(cur)*dens(cur))then
            if(cur.eq.1)then

                  flucount=flucount+1
                  if(ran.le.chance)then
                        cur=2    !1064nm
                  elseif(ran.gt.chance.and.ran.le..833)then
                        cur=3    !900nm
                  else
                        cur=4    !1300nm
                  end if
                        cost=2.*ran2(iseed)-1.
                        sint=(1.-cost*cost)
                        if(sint.le.0.)then
                              sint=0.
                        else
                              sint=sqrt(sint)
                        endif

                        phi=twopi*ran2(iseed)
                        sinp=sin(phi)
                        cosp=cos(phi)

                        nxp=sint*cosp
                        nyp=sint*sinp
                        nzp=cost
                        tauflag=.TRUE.
                        weight=1.
            end if
      end if
      

      
      end subroutine flurosub
