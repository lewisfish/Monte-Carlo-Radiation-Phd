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

   implicit none
   
   DOUBLE PRECISION :: Stratum_kappa, LiveEpi_kappa, PapDerm_kappa, RetDerm_kappa, HypoDerm_kappa
   DOUBLE PRECISION :: frac_H2O, nu_m, nu_Pd_Hb, nu_Rd_Hb
   DOUBLE PRECISION :: z, rho, mus_ref_500, mus_r, mus_m, a, f_ray, b, gam
   DOUBLE PRECISION :: nu_H2O ,nu_b, S, b_frac
   DOUBLE PRECISION :: mus1,mus2,mus3,mus4,mus5
   integer :: i,j,n

   hgg=0.7
   g2  = hgg**2.

   mus=0.
   nu_m = .001d0
   nu_Pd_Hb = 6.
   nu_Rd_hb = 4.5
   
   
!   open(12,file='abs.dat')
!do i=300,1000
!   wave=dble(i)
  
      
!      wave = HypoDerm_kappa+RetDerm_kappa+PapDerm_kappa+LiveEpi_kappa+Stratum_kappa
!write(12,*) i,HypoDerm_kappa, RetDerm_kappa, PapDerm_kappa, LiveEpi_kappa, Stratum_kappa, wave

!end do
!call exit(0)
!loop to set optical properties  
   do i=1,nzg
      z=zface(i)-zmax+zmax/nzg

      if(z.gt.zmax-0.02)then
         !Strat corenum
         rhokap(:,:,i,1) = Stratum_kappa
         albedo(:,:,i,1) = mus1/Stratum_kappa
      elseif(z.gt.zmax-0.1)then
         !Living Epidermis
         rhokap(:,:,i,1) = LiveEpi_kappa
         albedo(:,:,i,1) = mus2/LiveEpi_kappa
      elseif(z.gt.zmax-0.28)then
         !PaPillary Dermis
         rhokap(:,:,i,1) = PapDerm_kappa
         albedo(:,:,i,1) = mus3/PapDerm_kappa
      elseif(z.gt.zmax-2.1)then
         !Reticular Dermis
         rhokap(:,:,i,1) = RetDerm_kappa
         albedo(:,:,i,1) = mus4/RetDerm_kappa         
      elseif(z.lt.zmax-2.1)then
         !Hypodermis
         rhokap(:,:,i,1) = HypoDerm_kappa
         albedo(:,:,i,1) = mus5/HypoDerm_kappa
      end if
   end do

   end subroutine


   subroutine sample(array,size_of,cdf,wave,iseed)
!      
!  samples a random value from an array based upon its cdf     
!      
      implicit none
      
      integer, intent(IN)                :: iseed, size_of
      DOUBLE PRECISION,    intent(IN)    :: array(size_of,2), cdf(size_of)
      DOUBLE PRECISION,    intent(OUT)   :: wave
      
      real :: ran2
      DOUBLE PRECISION :: value
      integer :: nlow
      
      value = ran2(iseed)
      
      call search_1D(size(cdf),cdf,nlow,value)
      call lin_inter_1D(array,cdf,value,size(cdf),nlow,wave)
   
   end subroutine sample
   
   subroutine lin_inter_1D(array,cdf,value,length,nlow,y)
!
!  linear interpolates between values for an array and its cdf
!   
      implicit none
   
      DOUBLE PRECISION,    intent(OUT)  :: y
      integer, intent(IN)   :: length
      DOUBLE PRECISION,    intent(IN)   :: value, array(length, 2), cdf(length - 1)
      integer, intent(IN)   :: nlow
   
      y = array(nlow+1, 1) + (array(nlow+2, 1) - array(nlow + 1, 1)) * (value - cdf(nlow))/(cdf(nlow+1) - cdf(nlow))
   
   end subroutine lin_inter_1D
   
   subroutine lin_inter_2D(array,value,length,nlow,y)
!
!  linear interpolation for an array
!
      implicit none

      DOUBLE PRECISION,    intent(OUT)  :: y
      integer, intent(IN)   :: length
      DOUBLE PRECISION,    intent(IN)   :: value, array(length,2)
      integer, intent(IN)   :: nlow
   
      y = array(nlow,2) + (array(nlow+1,2) - array(nlow,2)) * (value - array(nlow,1))/(array(nlow+1,1) - array(nlow,1))
   
   end subroutine lin_inter_2D
   
   subroutine search_1D(length,array,nlow,value)
!
!  search by bisection for 1D array
!
      implicit none
      
      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      DOUBLE PRECISION,    intent(in)  :: array(length),value
      
      nup = length
      nlow = 1
      middle = int((nup+nlow)/2.)

      do while((nup - nlow).gt.1)
         middle = int((nup + nlow)/2.)
         if(value.gt.array(middle))then
            nlow = middle
         else
            nup = middle   
         end if
      end do
   end subroutine search_1D
   
   subroutine search_2D(length,array,nlow,value)
!
!  search by bisection for 2D array
!
      implicit none
      
      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      DOUBLE PRECISION,    intent(in)  :: array(length,2),value
      
      nup = length
      nlow = 1
      middle = int((nup+nlow)/2.)

      do while((nup - nlow).gt.1)
         middle = int((nup + nlow)/2.)
         if(value.gt.array(middle,1))then
            nlow = middle
         else
            nup = middle   
         end if
      end do
   end subroutine search_2D
   
   subroutine mk_cdf(array,cdf,length)
!
!  subroutine that creates cdf for an array of values.
!
      implicit none

      integer, intent(IN)    :: length
      DOUBLE PRECISION,    intent(IN)    :: array(length,2)
      DOUBLE PRECISION,    intent(INOUT) :: cdf(length)
      DOUBLE PRECISION                   :: summ
      integer                :: i,j
   
      do j=1,length-1
         summ=0.
         do i=1,j   
            summ=summ+0.5*(array(i+1,2)+array(i,2))*(array(i+1,1)-array(i,1))
         end do
         cdf(j)=summ      
      end do
      cdf=cdf/cdf(length-1)
   
   end subroutine mk_cdf

   DOUBLE PRECISION function water(wave)

   use iarray, only : water_array

      DOUBLE PRECISION ::  eps, wave
      integer :: nlow
            
      call search_2D(size(water_array,1),water_array,nlow,wave)
      call lin_inter_2D(water_array,wave,size(water_array,1),nlow,eps)

      water = eps !in cm-1

   end function water

   
   DOUBLE PRECISION function Oxy_Hb(conc, wave)
   
   use iarray, only : Oxy_Hb_array
   
      DOUBLE PRECISION :: conc, eps, wave
      integer :: nlow
         
      call search_2D(size(Oxy_Hb_array,1),Oxy_Hb_array,nlow,wave)
      call lin_inter_2D(Oxy_Hb_array,wave,size(Oxy_Hb_array,1),nlow,eps)
      
      Oxy_Hb = conc * ((eps)/66500.d0) !in cm-1

   end function Oxy_Hb

   DOUBLE PRECISION function Deoxy_Hb(conc, wave)

   use iarray, only : Deoxy_Hb_array

      DOUBLE PRECISION :: conc, eps, wave
      integer :: nlow
            
      call search_2D(size(Deoxy_Hb_array,1),Deoxy_Hb_array,nlow,wave)
      call lin_inter_2D(Deoxy_Hb_array,wave,size(Deoxy_Hb_array,1),nlow,eps)
      
      Deoxy_Hb = conc * ((eps)/66500.d0) !in cm-1

   end function Deoxy_Hb

   DOUBLE PRECISION function Carotene(conc, wave)

   use iarray, only : Carotene_array

      DOUBLE PRECISION :: conc, eps, wave
      integer :: nlow
      
      if(wave.gt.700)then
         eps=0.d0
      else
         call search_2D(size(Carotene_array,1),Carotene_array,nlow,wave)
         call lin_inter_2D(Carotene_array,wave,size(Carotene_array,1),nlow,eps)
         if(eps.lt.0)eps=0.
      end if
      Carotene = conc * ((eps)/537.d0) !in cm-1

   end function Carotene

   DOUBLE PRECISION function Bilirubin(conc, wave)

   use iarray, only : Bilirubin_array
   
      DOUBLE PRECISION :: conc, eps, wave
      integer :: nlow
           
           
      if(wave.gt.700)then
         eps=0.d0
      else
         call search_2D(size(Bilirubin_array,1),Bilirubin_array,nlow,wave)
         call lin_inter_2D(Bilirubin_array,wave,size(Bilirubin_array,1),nlow,eps)
         if(eps.lt.0)eps=0.
      end if
      
      Bilirubin = conc * ((eps)/585.d0) !in cm-1

   end function Bilirubin
   
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
   
end module ch_opt


!age 30

!vm=10 ; vHb=5 ; cbil=0.05 ; ccar=2.1e_4

![(.1-.3*10**(-4.)*wave) + .125*(wave/10.) * mua_base] * (1.-frac_H2O) + mua_H2O
