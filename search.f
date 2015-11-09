      subroutine search(value,array,length,Nlower)
      
      implicit none
      
      integer :: Nupper,Nlower,middle,length
      double precision :: value,array(length+1)
      
      Nupper=length
      Nlower=1
      middle=int((Nupper+Nlower)/2.)
      do while(Nupper-Nlower.gt.1.)
            middle=int((Nupper+Nlower)/2.)
            if(value.gt.array(middle))then
                  Nlower=middle 
            else
                  Nupper=middle
            end if

      end do
      end subroutine search
