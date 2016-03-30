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
   DOUBLE PRECISION    :: conc
   
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
   DOUBLE PRECISION :: z, mua, mus
   integer :: i,j
   
   !Strat Corneum sample
      !set mua
         call search_2D(size(mua_array,1),mua_array,nlow,wave)
         call lin_inter_2D(mua_array,wave,size(mua_array,1),nlow,mua)
      !set mus
         call search_2D(size(mus_array,1),mus_array,nlow,wave)
         call lin_inter_2D(mus_array,wave,size(mus_array,1),nlow,mus) 
   
   Stratum_kappa = mua + mus
   
   !Living Epidermis sample
      !set mua
         call search_2D(size(mua_array,1),mua_array,nlow,wave)
         call lin_inter_2D(mua_array,wave,size(mua_array,1),nlow,mua)
      !set mus
         call search_2D(size(mus_array,1),mus_array,nlow,wave)
         call lin_inter_2D(mus_array,wave,size(mus_array,1),nlow,mus) 
   
   LiveEpi_kappa = mua + mus
   
   !Pap Dermis sample
      !set mua
         call search_2D(size(mua_array,1),mua_array,nlow,wave)
         call lin_inter_2D(mua_array,wave,size(mua_array,1),nlow,mua)
      !set mus
         call search_2D(size(mus_array,1),mus_array,nlow,wave)
         call lin_inter_2D(mus_array,wave,size(mus_array,1),nlow,mus) 
   
   PapDerm_kappa = mua + mus
   
   !Ret Dermis Sample
      !set mua
         call search_2D(size(mua_array,1),mua_array,nlow,wave)
         call lin_inter_2D(mua_array,wave,size(mua_array,1),nlow,mua)
      !set mus
         call search_2D(size(mus_array,1),mus_array,nlow,wave)
         call lin_inter_2D(mus_array,wave,size(mus_array,1),nlow,mus) 
      
   RetDerm_kappa = mua + mus
   
   !Hypodermis smaple
      !set mua
         call search_2D(size(mua_array,1),mua_array,nlow,wave)
         call lin_inter_2D(mua_array,wave,size(mua_array,1),nlow,mua)
      !set mus
         call search_2D(size(mus_array,1),mus_array,nlow,wave)
         call lin_inter_2D(mus_array,wave,size(mus_array,1),nlow,mus) 
   
   HypoDerm_kappa = mua + mus

!loop to set optical properties  
   do i=1,nzg
      z=zface(i)-zmax+zmax/nzg
      
      if(z.gt.zmax-0.02)then
         !Strat corenum
         rhokap(:,:,i,1) =Stratum _kappa
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
end module ch_opt
