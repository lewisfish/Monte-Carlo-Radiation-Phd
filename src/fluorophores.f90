Module fluorophores

use subs

implicit none

CONTAINS

   DOUBLE PRECISION function NAD(wave)

   use iarray, only : NAD_array

      DOUBLE PRECISION ::  eps, wave, low
      integer :: nlow
      low = maxval(minval(nad_array,1))
      if(wave.lt.low)then
         nad=0.d0
      else
         call search_2D(size(NAD_array,1),NAD_array,nlow,wave)
         call lin_inter_2D(NAD_array,wave,size(NAD_array,1), nlow, NAD) !in cm-1
         if(nad .lt. 0.d0)nad=0.d0
      end if

   end function NAD
   
   DOUBLE PRECISION function NADH(wave)

   use iarray, only : NADH_array

      DOUBLE PRECISION ::  eps, wave
      integer :: nlow
            
      call search_2D(size(NADH_array,1),NADH_array,nlow,wave)
      call lin_inter_2D(NADH_array,wave,size(NADH_array,1), nlow, NADH) !in cm-1
      if(nadh .lt. 0.d0)nadh=0.d0
   end function NADH
   
   DOUBLE PRECISION function FAD(wave)

   use iarray, only : FAD_array

      DOUBLE PRECISION ::  eps, wave
      integer :: nlow
            
      call search_2D(size(FAD_array,1),FAD_array,nlow,wave)
      call lin_inter_2D(FAD_array,wave,size(FAD_array,1), nlow, FAD) !in cm-1
      if(fad .lt. 0.d0)fad=0.d0
   end function FAD
   
   DOUBLE PRECISION function RIBOFLAVIN(wave)

   use iarray, only : RIBOFLAVIN_array

      DOUBLE PRECISION ::  eps, wave
      integer :: nlow
            
      call search_2D(size(RIBOFLAVIN_array,1),RIBOFLAVIN_array,nlow,wave)
      call lin_inter_2D(RIBOFLAVIN_array,wave,size(RIBOFLAVIN_array,1), nlow, RIBOFLAVIN) !in cm-1
      if(RIBOFLAVIN .lt. 0.d0)RIBOFLAVIN=0.d0
   end function RIBOFLAVIN
            
   DOUBLE PRECISION function TYROSINE(wave)

   use iarray, only : TYROSINE_array

      DOUBLE PRECISION ::  eps, wave
      integer :: nlow
            
      call search_2D(size(TYROSINE_array,1),TYROSINE_array,nlow,wave)
      call lin_inter_2D(TYROSINE_array,wave,size(TYROSINE_array,1), nlow, TYROSINE) !in cm-1
      if(TYROSINE .lt. 0.d0)TYROSINE=0.d0
   end function TYROSINE
   
   DOUBLE PRECISION function tryptophan(wave)

   use iarray, only : tryptophan_array

      DOUBLE PRECISION ::  eps, wave
      integer :: nlow
            
      call search_2D(size(tryptophan_array,1),tryptophan_array,nlow,wave)
      call lin_inter_2D(tryptophan_array,wave,size(tryptophan_array,1), nlow, tryptophan) !in cm-1
      if(tryptophan .lt. 0.d0)tryptophan=0.d0
   end function tryptophan
end MODULE fluorophores
