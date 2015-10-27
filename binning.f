      subroutine binning(deposit,xcur,ycur,weight,ddx,fluro,
     +                  ddy,cbinsnum,zcur,ddz,cur)
      
      implicit none
      
      integer binz,binx,biny,cbinsnum,cur,rbins,binr
      real  deposit(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum),xcur,ycur
      real weight,zcur,ddx,ddy,ddz,r,ddr
      real fluro(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum,4)
      
!      if(cartbin.eq.1)then
      !cartesian bins
                  binx=floor(xcur/ddx)
                  biny=floor(ycur/ddy)

                  if(binx.lt.0.)binx=-1

                  if(biny.lt.0.)biny=-1

                  binz=floor(zcur/ddz)

                  if(binz.gt.cbinsnum)binz=-1
                  if(binz.lt.0.)binz=-1
            
                  if(cur.eq.1)then
            deposit(binx,biny,binz)=deposit(binx,biny,binz)+weight
                  else
            fluro(binx,biny,binz,cur)=fluro(binx,biny,binz,cur)+weight 
                  end if
!      else
      !cylindrical bins
      
!            binr=floor(r/ddr)
!            if(binr.lt.0.)binr=-1
!            if(binr.gt.rbins)binr=-1
!      
!      
!            binz=floor(zcur/ddz)
!            if(binz.gt.cbinsnum)binz=-1
!            if(binz.lt.0.)binz=-1
!      
!            if(cur.eq.1)then
!                  deposit(binr,binz)=deposit(binr,binz)+weight
!            else
!                  fluro(binr,binz)=fluro(binr,binz)+weight
!            end if

      end subroutine binning
