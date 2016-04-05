module scatt_method


implicit none
save

contains
   
   subroutine russian_roul(nscatt, iseed, tflag, weight)
      
      use opt_prop
      use binning_mod
      use stokes_mod
      use constants, only : TERMINATE, CHANCE
      
      implicit none
      
      integer :: iseed
      DOUBLE PRECISION :: nscatt
      logical :: tflag
      DOUBLE PRECISION :: ran, weight, absorb
      real :: ran2
      
      ran = ran2(iseed)
      
!      weight = weight * albedo
      absorb = weight * (mua)/kappa
      call binning(absorb)
      nscatt = nscatt + 1
      
      call stokes(iseed)
      
      if(weight .lt. TERMINATE)then
         if(ran .le. CHANCE)then
            weight = weight/CHANCE
         else
            absorb = weight
            call binning(absorb)
            weight = 0.
            tflag = .TRUE.
         end if
      end if
   
   end subroutine russian_roul
   
   subroutine normal_S_A(nscatt, iseed, tflag)

!      use opt_prop, only : albedo
      use constants, only : acount
      use binning_mod      
      use stokes_mod      
      implicit none
      
      DOUBLE PRECISION :: ran
      real :: ran2
      real*8 :: nscatt
      integer :: iseed
      logical :: tflag
      
      ran = ran2(iseed)
      
!!      if(ran .lt. albedo)then
!         call stokes(iseed)
!         nscatt = nscatt + 1
!      else
!         tflag = .TRUE.
!         call binning
!         acount = acount + 1
!      end if
      
   end subroutine normal_S_A
   
   subroutine normal_S_A_F(nscatt, iseed, hgghold, tflag)

      use opt_prop    
      use stokes_mod
      use binning_mod
      use iarray, only : e_cdf_c, excite_array_c
      use constants, only : acount, fcount
      use ch_opt, only : init_opt, sample
      
      implicit none
      
      DOUBLE PRECISION :: nscatt, ran, hgghold
      real :: ran2
      integer :: iseed
      logical :: tflag
      
      ran = ran2(iseed)
      
!!      if(ran .lt. albedo)then
!         call stokes(iseed)
!         nscatt = nscatt + 1
!      elseif(ran .lt. (muaf+mus)/kappa)then
!         call sample(excite_array_c,size(e_cdf_c),e_cdf_c,wave,iseed)
!         call init_opt
!         fcount = fcount+1
!         hgg = 0.
!         g2 = 0.
!         call stokes(iseed)
!         hgg = hgghold
!         g2 = hgg**2.
!      else
!         tflag = .TRUE.
!         call binning
!         acount = acount + 1
!      end if
      
   end subroutine normal_S_A_F
end module scatt_method 
