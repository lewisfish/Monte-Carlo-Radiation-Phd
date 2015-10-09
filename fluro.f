      subroutine flurosub(flucount,rhokap,iseed,xcell,ycell,zcell,cur
     +                   ,kappa2,nxp,nyp,sint,cost,sinp,phi,hgg,g2,pi,
                         twopi)
      
      implicit none
      
      include 'grid.txt'
    !subroutine to check if photon is in fluro material, then release a fluro photon if appropriate
      integer xcell,ycell,zcell,cur,iseed,flucount
      real ran2,kappa2,chance
      
      chance=0.05
      
      if(rhokap(xcell,ycell,zcell,cur).eq.kappa2)then
            if(cur.eq.1)then
                  if(ran2(iseed).lt.chance)then
                        flucount=flucount+1
                        cur=2
                        fflag=1
                    call stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
     +                  hgg,g2,pi,twopi,iseed,cur,fflag)
                  end if
            end if
      end if
      

      
      end subroutine flurosub
