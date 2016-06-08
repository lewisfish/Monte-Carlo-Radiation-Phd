Module fluorophores

use subs

implicit none

CONTAINS
   
   DOUBLE PRECISION function NADH(wave)

   use iarray, only : NADH_array

      DOUBLE PRECISION ::  eps, wave, low, high 
      integer :: nlow, sizeof
      
      sizeof = size(nadh_array,1)
      low    = nadh_array(1,1)
      high   = nadh_array(sizeof,1)
      
      if(wave .le. low .or. wave .gt. high)then
         nadh = 0.0d0
      else
         call search_2D(size(NADH_array,1),NADH_array,nlow,wave)
         call lin_inter_2D(NADH_array,wave,size(NADH_array,1), nlow, NADH) !in cm-1
         if(nadh .lt. 0.d0)nadh=0.d0
      end if
   end function NADH
   
      
   DOUBLE PRECISION function RIBOFLAVIN(wave)

   use iarray, only : RIBOFLAVIN_array

      DOUBLE PRECISION ::  eps, wave, low, high 
      integer :: nlow, sizeof
      
      sizeof = size(riboflavin_array,1)
      low    = riboflavin_array(1,1)
      high   = riboflavin_array(sizeof,1)
      
      if(wave .le. low .or. wave .gt. high)then
         riboflavin = 0.0d0
      else
         call search_2D(size(RIBOFLAVIN_array,1),RIBOFLAVIN_array,nlow,wave)
         call lin_inter_2D(RIBOFLAVIN_array,wave,size(RIBOFLAVIN_array,1), nlow, RIBOFLAVIN) !in cm-1
         if(RIBOFLAVIN .lt. 0.d0)RIBOFLAVIN=0.d0
      end if
   end function RIBOFLAVIN
            
   DOUBLE PRECISION function TYROSINE(wave)

      use iarray, only : TYROSINE_array

      DOUBLE PRECISION ::  eps, wave, low, high 
      integer :: nlow, sizeof
      
      sizeof = size(TYROSINE_array,1)
      low    = TYROSINE_array(1,1)
      high   = TYROSINE_array(sizeof,1)
      
      if(wave .le. low .or. wave .gt. high)then
         tyrosine = 0.0d0
      else
         call search_2D(size(TYROSINE_array,1),TYROSINE_array,nlow,wave)
         call lin_inter_2D(TYROSINE_array,wave,size(TYROSINE_array,1), nlow, TYROSINE) !in cm-1
         if(TYROSINE .lt. 0.d0)TYROSINE=0.d0
      end if
   end function TYROSINE
end MODULE fluorophores
