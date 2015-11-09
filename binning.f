      subroutine binning(deposit,xcur,ycur,weight,ddx,fluro,
     +                  ddy,cbinsnum,zcur,ddz,cur)
      
      implicit none
      
      integer binz,binx,biny,cbinsnum,cur
      double precision  deposit(-1:cbinsnum,-1:cbinsnum,
     +      -1:cbinsnum)
      double precision weight,zcur,ddx,ddy,ddz,xcur,ycur
      double precision fluro(-1:cbinsnum,-1:cbinsnum,-1:
     +      cbinsnum,4)

                  binx=floor(xcur/ddx)
                  biny=floor(ycur/ddy)

                  if(binx.lt.0.)binx=-1

                  if(biny.lt.0.)biny=-1

                  binz=floor(zcur/ddz)
!                  if(cur.eq.4)then
!                  if(binz.gt.170)print*,cbinsnum,binz
!                  end if
                  if(binz.gt.cbinsnum)binz=-1
                  if(binz.lt.0.)binz=-1

                  if(cur.eq.1)then
            deposit(binx,biny,binz)=deposit(binx,biny,binz)+weight
                  else
            fluro(binx,biny,binz,cur)=fluro(binx,biny,binz,cur)+weight
                  end if


      end subroutine binning
