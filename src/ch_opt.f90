MODULE ch_opt

implicit none

CONTAINS
   
   subroutine init_opt
!
!  subroutine to set optical properties
!

   use iarray, only : mus_array,mua_array
   use opt_prop
   
   implicit none
   
   integer :: nlow
   DOUBLE PRECISION :: conc
   
   !set mua_skin
   call search_2D(size(mua_array,1),mua_array,nlow,wave)
   call lin_inter_2D(mua_array,wave,size(mua_array,1),nlow,mua)
   !set mus_skin
   call search_2D(size(mus_array,1),mus_array,nlow,wave)
   call lin_inter_2D(mus_array,wave,size(mus_array,1),nlow,mus)   

   hgg=0.7
   g2  = hgg**2.
   kappa  = mus + mua
   albedo = (mus)/kappa

   end subroutine init_opt


   subroutine opt_set(zmax, wave)
!
!
!
   use constants, only : nzg,nyg,nxg
   use iarray,    only : zface, rhokap

   implicit none
   
   DOUBLE PRECISION, intent(IN) :: zmax, wave
   DOUBLE PRECISION :: Stratum_kappa, LiveEpi_kappa, PapDerm_kappa, RetDerm_kappa, HypoDerm_kappa
   DOUBLE PRECISION :: frac_H2O, nu_m, nu_Pd_Hb, nu_Rd_Hb, conc, car_conc
   DOUBLE PRECISION :: z, mua, mus
   integer :: i,j
   
   !Strat Corneum sample
      !set mua
         mua = ((0.1 - 0.3*10**(-4.) * wave) + 0.125 * (wave/10.) * Base(wave)) * (1. - frac_H2O) + water(wave)   
      !set mus
         
   Stratum_kappa = mua + mus
   
   !Living Epidermis sample
      !set mua
         mua = (nu_m * (Eumel(wave) + Pheomel(wave)) + (1. + nu_m) * Carotene(conc, wave)) * (1.-frac_H2O) + water(wave) 
      !set mus
         
   LiveEpi_kappa = mua + mus
   
   !Pap Dermis sample
      !set mua
         mua = (nu_Pd_Hb * (Oxy_Hb(conc, wave) + Deoxy_Hb(conc, wave) + Bilirubin(conc, wave) + &
               Carotene(conc, wave) + (1.- Car_conc) * base(wave))) * (1.-frac_H2O) + water(wave) 
      !set mus
         
   PapDerm_kappa = mua + mus
   
   !Ret Dermis Sample
      !set mua
         mua = (nu_Rd_Hb * (Oxy_Hb(conc, wave) + Deoxy_Hb(conc, wave) + Bilirubin(conc, wave) + &
               Carotene(conc, wave) + (1.- Car_conc) * base(wave))) * (1.-frac_H2O) + water(wave)
      !set mus
      
   RetDerm_kappa = mua + mus
   
   !Hypodermis smaple
      !set mua
         mua = water(wave)
      !set mus
   
   HypoDerm_kappa = mua + mus

!loop to set optical properties  
   do i=1,nzg
      z=zface(i)-zmax+zmax/nzg
      
      if(z.gt.zmax-0.02)then
         !Strat corenum
         rhokap(:,:,i,1) =Stratum_kappa
      elseif(z.gt.zmax-0.1)then
         !Living Epidermis
         rhokap(:,:,i,1) = LiveEpi_kappa
      elseif(z.gt.zmax-0.28)then
         !PaPillary Dermis
         rhokap(:,:,i,1) = PapDerm_kappa
      elseif(z.gt.zmax-2.1)then
         !Reticular Dermis
         rhokap(:,:,i,1) = RetDerm_kappa
      elseif(z.lt.zmax-2.1)then
         !Hypodermis
         rhokap(:,:,i,1) = HypoDerm_kappa
      end if
   end do

!   INQUIRE(iolength = i) array
!   print*,'irec ',i
!   open(62,file='arraytest.dat',access='direct',form='unformatted',recl=i)
!   write(62,rec=1) array
!   close(62)
   
!   open(12,file='test.dat')
!   do i=1,nzg
!      write(12,*) (array(i,50,j),j=1,nzg)
!   end do
   
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

   use iarray, only : water_cdf, water_array

      DOUBLE PRECISION ::  eps, wave
      integer :: nlow
            
      call search_2D(size(water_cdf),water_array,nlow,wave)
      call lin_inter_2D(water_array,wave,size(water_cdf),nlow,eps)
      
      water = eps

   end function water

   
   DOUBLE PRECISION function Oxy_Hb(conc, wave)
   
   use iarray, only : Oxy_Hb_cdf, Oxy_Hb_array
   
      DOUBLE PRECISION :: conc, eps, wave
      integer :: nlow
         
      call search_2D(size(Oxy_Hb_cdf),Oxy_Hb_array,nlow,wave)
      call lin_inter_2D(Oxy_Hb_array,wave,size(Oxy_Hb_cdf),nlow,eps)
      
      Oxy_Hb = conc * ((eps)/66500.d0)

   end function Oxy_Hb

   DOUBLE PRECISION function Deoxy_Hb(conc, wave)

   use iarray, only : Deoxy_Hb_cdf, Deoxy_Hb_array

      DOUBLE PRECISION :: conc, eps, wave
      integer :: nlow
            
      call search_2D(size(Deoxy_Hb_cdf),Deoxy_Hb_array,nlow,wave)
      call lin_inter_2D(Deoxy_Hb_array,wave,size(Deoxy_Hb_cdf),nlow,eps)
      
      Deoxy_Hb = conc * ((eps)/66500.d0)

   end function Deoxy_Hb

   DOUBLE PRECISION function Carotene(conc, wave)

   use iarray, only : Carotene_cdf, Carotene_array

      DOUBLE PRECISION :: conc, eps, wave
      integer :: nlow
            
      call search_2D(size(Carotene_cdf),Carotene_array,nlow,wave)
      call lin_inter_2D(Carotene_array,wave,size(Carotene_cdf),nlow,eps)
      
      Carotene = conc * ((eps)/537.d0)

   end function Carotene

   DOUBLE PRECISION function Bilirubin(conc, wave)

   use iarray, only : Bilirubin_cdf, Bilirubin_array
   
      DOUBLE PRECISION :: conc, eps, wave
      integer :: nlow
      
      call search_2D(size(Bilirubin_cdf),Bilirubin_array,nlow,wave)
      call lin_inter_2D(Bilirubin_array,wave,size(Bilirubin_cdf),nlow,eps)
      
      Bilirubin = conc * ((eps)/585.d0)

   end function Bilirubin
   
   DOUBLE PRECISION function Base(wave)

      DOUBLE PRECISION :: wave

      Base = 1.78*10.**7. * wave**(-3.255)
   
   end function Base
   
   DOUBLE PRECISION function Eumel(wave)

      DOUBLE PRECISION :: wave

      Eumel = 6.6*10.**10. * wave**(-3.33)
   
   end function Eumel
   
   DOUBLE PRECISION function Pheomel(wave)

      DOUBLE PRECISION :: wave

      Pheomel = 2.9*10.**14. * wave**(-4.75)
   
   end function Pheomel
end module ch_opt


!age 30

!vm=10 ; vHb=5 ; cbil=0.05 ; ccar=2.1e_4

![(.1-.3*10**(-4.)*wave) + .125*(wave/10.) * mua_base] * (1.-frac_H2O) + mua_H2O
