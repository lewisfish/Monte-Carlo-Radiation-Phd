      subroutine flurosub(flucount,rhokap,iseed,xcell,ycell,zcell,cur
     +                   ,kappa2)
      
      implicit none
      
      include 'grid.txt'
    
      real ran2,kappa2,chance
      integer xcell,ycell,zcell,cur,iseed,flucount
      
      chance=0.05
      
      if(rhokap(xcell,ycell,zcell,cur).eq.kappa2)then
            if(cur.eq.1)then
                  if(ran2(iseed).lt.chance)then
                        flucount=flucount+1
                        cur=2
                  end if
            end if
      end if
      
      end subroutine flurosub
