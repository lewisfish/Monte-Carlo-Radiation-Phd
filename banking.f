      subroutine banking(bank,xcur,ycur,zcur,weight,nxp,nyp,nzp
     +                  ,WH,WL,totcount,loopflag,ibank)
      
      implicit none
      
      real :: bank(0:10000,1:7,0:1),WH,WL
      real :: xcur,ycur,zcur,weight,nxp,nyp,nzp
      integer :: counter,i,totcount,ibank
      logical :: loopflag
      
      counter=0
      if(weight.gt.WH)then
      !split particles
            do while(weight.gt.WH)
                  weight=weight/2.
                  counter=counter+2
            end do
            do i=totcount,totcount+counter
            
                  bank(totcount+counter,1,ibank)=xcur
                  bank(totcount+counter,2,ibank)=ycur
                  bank(totcount+counter,3,ibank)=zcur
                  bank(totcount+counter,4,ibank)=nxp
                  bank(totcount+counter,5,ibank)=nyp
                  bank(totcount+counter,6,ibank)=nzp
                  bank(totcount+counter,7,ibank)=weight
            
            end do
            totcount=totcount+counter
            counter=0
      elseif(weight.lt.WL)then
      ! play russian roulette
      end if
      

      
      
      end subroutine banking
