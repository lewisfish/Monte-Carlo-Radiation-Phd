MODULE ch_opt

implicit none

CONTAINS
   
   subroutine init_opt
   
   use iarray, only : mua_array,mus_array,excite_array,fluro_array
   use opt_prop
   
   implicit none
   
   integer :: nlow

   !set mua_skin
   call search_2D(size(fluro_array,1),fluro_array,nlow,wave)
   call lin_inter_2D(fluro_array,wave,size(fluro_array,1),nlow,mua)
   if(mua.lt.0.)mua=0.
   !set mua_fluro
!   call search_2D(size(excite_array,1),excite_array,nlow,wave)
!   call lin_inter_2D(excite_array,wave,size(excite_array,1),nlow,muaf)   

!   set mus
   call search_2D(size(mus_array,1),mus_array,nlow,wave)
   call lin_inter_2D(mus_array,wave,size(mus_array,1),nlow,mus)
!   set g and hgg
!   hgg = 0.62 + 0.29 * 10.**(-3.) * wave
   hgg=0.7
   g2  = hgg**2.
!   mua = 80.
!   mus = 250.
!   muai=mus/999.
   kappa  = mus + mua + (mus/999.) 
   albedo = mus/kappa
!   print*,mus,mua,albedo
!   mus / kappa

   end subroutine init_opt

   subroutine sample(array,size_of,cdf,wave,iseed)
      
      integer, intent(IN)    :: iseed,size_of
      real,    intent(IN)    :: array(size_of,2),cdf(size_of)
      real,    intent(INOUT) :: wave

      real :: ran2,value
      integer :: nlow
      
      value = ran2(iseed)
      
      call search_1D(size(cdf),cdf,nlow,value)
      call lin_inter_1D(array,cdf,value,size(cdf),nlow,wave)
   
   end subroutine sample
   
   subroutine lin_inter_1D(array,cdf,value,length,nlow,y)
   
      real, intent(OUT)  :: y
      integer, intent(IN) :: length
      real, intent(IN)   :: value,array(length,2),cdf(length-1)
      integer,intent(IN) :: nlow
   
      y = array(nlow+1,1) + (array(nlow+2,1) - array(nlow+1,1)) * (value - cdf(nlow))/(cdf(nlow+1) - cdf(nlow))
   
   end subroutine lin_inter_1D
   
   subroutine lin_inter_2D(array,value,length,nlow,y)
   
      real, intent(OUT)  :: y
      integer, intent(IN) :: length
      real, intent(IN)   :: value,array(length,2)
      integer,intent(IN) :: nlow
   
      y = array(nlow,2) + (array(nlow+1,2) - array(nlow,2)) * (value - array(nlow,1))/(array(nlow+1,1) - array(nlow,1))
   
   end subroutine lin_inter_2D
   
   subroutine search_1D(length,array,nlow,value)
      
      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      real, intent(in)     :: array(length),value
      
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
      
      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      real, intent(in)     :: array(length,2),value
      
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

   integer, intent(IN) :: length
   real, intent(IN)    :: array(length,2)
   real, intent(INOUT) :: cdf(length)
   real                :: summ
   integer             :: i,j
   
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
