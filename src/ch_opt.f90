MODULE ch_opt

implicit none

CONTAINS
   
   subroutine init_opt(conc)
!
!  subroutine to set optical properties
!

   use iarray, only : mus_array, mua_array!, conc
   use opt_prop
   use fluorophores
   use subs
   
   implicit none
   
   integer :: nlow
   double precision :: conc(3)
   
   !set mua_skin
   call search_2D(size(mua_array,1),mua_array,nlow,wave)
   call lin_inter_2D(mua_array,wave,size(mua_array,1),nlow,mua)
   
   !set mus_skin
   call search_2D(size(mus_array,1),mus_array,nlow,wave)
   call lin_inter_2D(mus_array,wave,size(mus_array,1),nlow,mus)   
   
   !set mua_nadh, riboflavin, and tyrosine
   mua_nadh   = 2.3d0 * (nadh(wave) * conc(1))
   mua_ribo   = 2.3d0 * (riboflavin(wave) * conc(2))
   mua_tyro   = 2.3d0 * (tyrosine(wave) * conc(3))
   
   hgg    = 0.7d0
   g2     = hgg**2.
   kappa  = mus + mua + mua_nadh + mua_ribo + mua_tyro
   albedo = (mus) / kappa

   end subroutine init_opt
end module ch_opt
