module scatt_method

save
implicit none

contains
   
   subroutine russian_roul(nscatt, iseed, weight, tflag)
      
      use opt_prop
      use binning_mod
      
      implicit none
      
      integer :: iseed
      real*8 :: nscatt
      logical :: tflag
      real :: ran2, weight
   
      ran = ran2(iseed)
      
      weight = weight * albedo
      absorb = weight * (mua)/kappa
      !call binning()
      nscatt = nscatt + 1
      
      call stokes(iseed)
      
      if(weight .lt. terminate)then
         if(ran .le. chance)then
            weight = weight/chance
         else
            absorb = weight
            !call binning()
            weight = 0.
            tflag = .TRUE.
         end if
      end if
   
   end subroutine russian_roul
   
   subroutine normal_S_A(nscatt, iseed, counter2, tflag)

      use opt_prop, only : albedo      
      use stokes_mod      
      implicit none
      
      real :: ran2,ran
      real*8 :: nscatt
      integer :: iseed, counter2
      
      ran = ran2(iseed)
      
      if(ran .lt. albedo)then
         call stokes(iseed)
         nscatt = nscatt + 1
      else
         tflag = .TRUE.
         counter2 = counter2 + 1
      end if
      
   end subroutine normal_S_A
   
   subroutine normal_S_A_F(nscatt, iseed, counter2, counter1, hgghold, tflag)

      use opt_prop    
      use stokes_mod
      use ch_opt, only : init_opt, sample
      
      implicit none
      
      real*8 :: nscatt
      real :: ran2, ran, hgghold
      integer :: iseed, counter1, counter2
      logical :: tflag
      
      ran = ran2(iseed)
      
      if(ran .lt. albedo)then
         call stokes(iseed)
         nscatt = nscatt + 1
      elseif(ran .lt. (muaf+mus)/kappa)then
         call sample(excite_array,size(e_cdf),e_cdf,wave,iseed)
         call init_opt
!         fluro_pos(xcell,ycell,zcell)=fluro_pos(xcell,ycell,zcell)+1
         counter2=counter2+1
         hgg=0.
         g2=0.
         call stokes(iseed)
         hgg=hgghold
         g2=hgg**2.
         counter2=counter2+1
      else
         tflag = .TRUE.
         counter1=counter1+1
      end if
      
   end subroutine normal_S_A
end module scatt_method 
