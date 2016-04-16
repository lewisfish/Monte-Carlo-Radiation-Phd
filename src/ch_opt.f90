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

   implicit none
   
   DOUBLE PRECISION :: z
   integer :: i,j
   wave=500.d0
!loop to set optical properties in z dir. uniform in x,y 
   do i=1,nzg
      z=zface(i)-zmax+zmax/nzg

      if(z.gt.zmax-0.02)then
         !Strat corenum
         rhokap(:,:,i,1) = Stratum(wave)
!         albedo(:,:,i,1) = mus1/Stratum_kappa
      elseif(z.gt.zmax-0.1)then
         !Living Epidermis
         rhokap(:,:,i,1) = Epidermis(wave)
!         albedo(:,:,i,1) = mus2/LiveEpi_kappa
      elseif(z.gt.zmax-0.28)then
         !PaPillary Dermis
         rhokap(:,:,i,1) = Pap_Dermis(wave)
!         albedo(:,:,i,1) = mus3/PapDerm_kappa
      elseif(z.gt.zmax-2.1)then
         !Reticular Dermis
         rhokap(:,:,i,1) = Ret_Dermis(wave)
!         albedo(:,:,i,1) = mus4/RetDerm_kappa         
      elseif(z.lt.zmax-2.1)then
         !Hypodermis
         rhokap(:,:,i,1) = Hypo_Dermis(wave)
!         albedo(:,:,i,1) = mus5/HypoDerm_kappa
      end if
   end do

   inquire(iolength=i)rhokap(:,:,:,1)
   open(62,file='skintest.dat')
   do i=1,nxg
   write(62,*)(rhokap(i,100,j,1),j=1,nzg)
   end do
   close(62)
   call exit(0)
   end subroutine
end module ch_opt
