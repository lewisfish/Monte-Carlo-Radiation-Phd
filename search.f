      subroutine search(value,array,length,iseed,Nlower)
      
      implicit none
      
      real :: value,array(1:37),ran2
      integer :: Nupper,Nlower,middle,length,iseed
      
      
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
