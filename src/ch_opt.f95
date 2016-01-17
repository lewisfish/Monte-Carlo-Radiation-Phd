MODULE ch_opt

implicit none
save

CONTAINS
   subroutine init_opt
   
   use iarray, only : mua_array,mus_array
   use opt_prop
   
   implicit none
   
   integer :: Nlower
   real    :: wave
   
   wave=500.
   
   !set mua
   call search(wave,mua_array,size(mua_array,1),Nlower)
   call inter_lin(mua_array,size(mua_array,1),Nlower,wave,mua)
   
!   set mus
   call search(wave,mus_array,size(mus_array,1),Nlower)
   call inter_lin(mus_array,size(mus_array,1),Nlower,wave,mus)

!   set g and hgg
   hgg = 0.62 + 0.29 * 10.**(-3.) * wave
   g2  = hgg**2
   
   kappa  = mua + mus 
   albedo = mus / kappa
   
   end subroutine init_opt
   
   subroutine search(value,array,length,Nlower)

      implicit none

      integer          :: Nupper,Nlower,middle,length,j
      real, intent(IN) :: value,array(length,2)

      Nupper=length
      Nlower=1
      middle=int((Nupper+Nlower)/2.)

      do while(Nupper-Nlower.gt.1.)
         middle=int((Nupper+Nlower)/2.)
         if(value.gt.array(middle,1))then
            Nlower=middle 
         else
            Nupper=middle
         end if

      end do

      end subroutine search
      
      subroutine inter_lin(array,length,i,value,y)
      
      implicit none
      
      integer, intent(IN)  :: length,i
      real,    intent(IN)  :: array(length,2),value
      real,    intent(OUT) :: y
      
      y=array(i,2)+(array(i+1,2)-array(i,2))*((value-array(i,1))/(array(i+1,1)-array(i,1)))

      end subroutine inter_lin
      
end MODULE ch_opt


