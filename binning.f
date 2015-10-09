      subroutine binning(deposit,xcur,ycur,weight,ddx,fluro,
     +                  ddy,cbinsnum,zcur,ddz,cur)
      
      implicit none
      
      integer binz,binx,biny,cbinsnum,cur
      real  deposit(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum),xcur,ycur
      real weight,zcur,ddx,ddy,ddz
      real fluro(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum)
      
      
            binx=floor(xcur/ddx)
            biny=floor(ycur/ddy)
!            if(binx.gt.cbinsnum)binx=-1
            if(binx.lt.0.)binx=-1
!            if(biny.gt.cbinsnum)biny=-1
            if(biny.lt.0.)biny=-1

            binz=floor(zcur/ddz)

            if(binz.gt.cbinsnum)binz=-1
            if(binz.lt.0.)binz=-1
            
            if(cur.eq.1)then
            deposit(binx,biny,binz)=deposit(binx,biny,binz)+weight
            else
            fluro(binx,biny,binz)=fluro(binx,biny,binz)+weight 
            end if

      end subroutine binning
