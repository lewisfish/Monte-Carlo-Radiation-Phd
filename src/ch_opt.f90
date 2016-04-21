MODULE ch_opt

implicit none

CONTAINS
   
   subroutine init_opt
!
!  subroutine to set optical properties
!
   use opt_prop
   
   implicit none
   
   integer :: nlow
   DOUBLE PRECISION :: conc
   
   !set mua_skin
!   call search_2D(size(mua_array,1),mua_array,nlow,wave)
!   call lin_inter_2D(mua_array,wave,size(mua_array,1),nlow,mua)
!   !set mus_skin
!   call search_2D(size(mus_array,1),mus_array,nlow,wave)
!   call lin_inter_2D(mus_array,wave,size(mus_array,1),nlow,mus)   

!   hgg=0.7
!   g2  = hgg**2.
!   kappa  = mus + mua
!   albedo = (mus)/kappa

   end subroutine init_opt


   subroutine opt_set()
!
!
!
   use constants, only : nzg,nyg,nxg,zmax
   use iarray,    only : zface, rhokap, albedo
   use opt_prop
   use skin_prop
   use photon_vars, only : zp

   implicit none
   
   DOUBLE PRECISION :: z, kappa_tmp
   integer :: i,j
!   wave=500.d0
!loop to set optical properties in z dir. uniform in x,y 
!   do i=1,nzg
!      z=zface(i)-zmax+zmax/nzg

!      if(z.gt.zmax-0.02)then
         !Strat corenum                                  ! 0 means mua + mus  |  1 means mus
         kappa_tmp = Stratum(wave, 200, 0)
         rhokap(:,:,200,1) = kappa_tmp
         albedo(:,:,200,1) = Stratum(wave, 200, 1) / kappa_tmp
!      elseif(z.gt.zmax-0.1)then
         !Living Epidermis
         kappa_tmp = Epidermis(wave, 198, 0)
         rhokap(:,:,197:199,1) = kappa_tmp
         albedo(:,:,197:199,1) = Epidermis(wave, 198, 1) / kappa_tmp
!      elseif(z.gt.zmax-0.28)then
         if(zp .lt. zmax-0.05)then
         !PaPillary Dermis
         kappa_tmp = Pap_Dermis(wave, 194, 0)
         rhokap(:,:,190:196,1) = kappa_tmp
         albedo(:,:,190:196,1) = Pap_Dermis(wave, 194, 1) / kappa_tmp
!      elseif(z.gt.zmax-2.1)then
         elseif(zp .lt. zmax-.1d0)then
         !Reticular Dermis
         kappa_tmp = Ret_Dermis(wave, 140, 0)
         rhokap(:,:,117:189,1) = kappa_tmp
         albedo(:,:,117:189,1) = Ret_Dermis(wave, 140, 1) / kappa_tmp         
!      elseif(z.lt.zmax-2.1)then
         elseif(zp .lt. zmax-1.0d0)then
         !Hypodermis
         kappa_tmp = Hypo_Dermis(wave, 50, 0)
         rhokap(:,:,1:116,1) = kappa_tmp
         albedo(:,:,1:116,1) = Hypo_Dermis(wave, 50, 1) / kappa_tmp 
      end if
!   end do

   end subroutine opt_set
end module ch_opt
