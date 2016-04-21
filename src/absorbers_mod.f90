Module absorbers

use subs

implicit none

CONTAINS

   DOUBLE PRECISION function water(wave)

   use iarray, only : water_array

      DOUBLE PRECISION ::  eps, wave
      integer :: nlow
            
      call search_2D(size(water_array,1),water_array,nlow,wave)
      call lin_inter_2D(water_array,wave,size(water_array,1),nlow,eps)

      water = eps !in cm-1

   end function water

   DOUBLE PRECISION function Oxy_Hb(wave)
   
   use iarray, only : Oxy_Hb_array
   
      DOUBLE PRECISION :: eps, wave
      integer :: nlow
         
      call search_2D(size(Oxy_Hb_array,1),Oxy_Hb_array,nlow,wave)
      call lin_inter_2D(Oxy_Hb_array,wave,size(Oxy_Hb_array,1),nlow,eps)
      
      Oxy_Hb = ((eps)/66500.d0) !in cm-1

   end function Oxy_Hb

   DOUBLE PRECISION function Deoxy_Hb(wave)

   use iarray, only : Deoxy_Hb_array

      DOUBLE PRECISION :: eps, wave
      integer :: nlow
            
      call search_2D(size(Deoxy_Hb_array,1),Deoxy_Hb_array,nlow,wave)
      call lin_inter_2D(Deoxy_Hb_array,wave,size(Deoxy_Hb_array,1),nlow,eps)
      
      Deoxy_Hb = ((eps)/66500.d0) !in cm-1

   end function Deoxy_Hb

   DOUBLE PRECISION function Carotene(wave)

   use iarray, only : Carotene_array

      DOUBLE PRECISION :: eps, wave
      integer :: nlow
      
      if(wave.gt.700)then
         eps=0.d0
      else
         call search_2D(size(Carotene_array,1),Carotene_array,nlow,wave)
         call lin_inter_2D(Carotene_array,wave,size(Carotene_array,1),nlow,eps)
         if(eps.lt.0)eps=0.
      end if
      Carotene = ((eps)/537.d0) !in cm-1

   end function Carotene

   DOUBLE PRECISION function Bilirubin(wave)

   use iarray, only : Bilirubin_array
   
      DOUBLE PRECISION :: eps, wave
      integer :: nlow
           
           
      if(wave.gt.700)then
         eps=0.d0
      else
         call search_2D(size(Bilirubin_array,1),Bilirubin_array,nlow,wave)
         call lin_inter_2D(Bilirubin_array,wave,size(Bilirubin_array,1),nlow,eps)
         if(eps.lt.0)eps=0.
      end if
      
      Bilirubin = ((eps)/585.d0) !in cm-1

   end function Bilirubin
   
   DOUBLE PRECISION function fat(wave)
   
   use iarray, only : fat_array
   
   DOUBLE PRECISION, intent(in) :: wave
   integer                      :: nlow
   
   if(wave .le. 429.d0)then
      fat = 0.0d0
   else
         call search_2D(size(fat_array,1),fat_array,nlow,wave)
         call lin_inter_2D(fat_array,wave,size(fat_array,1),nlow,fat)
      if(fat .lt. 0.d0) fat =0.d0
   end if
   
   end function fat   
   
   DOUBLE PRECISION function Base(wave)

      DOUBLE PRECISION :: wave

      Base = 7.84*10.**7. * wave**(-3.255) !in cm-1

   end function Base
   
   DOUBLE PRECISION function Eumel(wave)

      DOUBLE PRECISION :: wave

      Eumel = 6.6*10.**11. * wave**(-3.33) !in cm-1
   
   end function Eumel
   
   DOUBLE PRECISION function Pheomel(wave)

      DOUBLE PRECISION :: wave

      Pheomel = 2.9*10.**15. * wave**(-4.75) !in cm-1
   
   end function Pheomel
end MODULE absorbers
