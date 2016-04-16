MODULE skin_prop

use opt_prop, only : hgg
use absorbers

implicit none

CONTAINS
   DOUBLE PRECISION function stratum(wave)
   
      DOUBLE PRECISION, intent(in) :: wave
      DOUBLE PRECISION             :: mua, rho, gam, mus_ref_500, mus_r, mus_m, mus, nu_H2O


      if(wave .ge. 0.0d0)then
      !set mua
      ! formula from MM02
      nu_H2O = .5
         mua = nu_H2O * water(wave) + (1.0d0 - nu_H2O) * base(wave)    !in cm-1   
      !set mus
      ! formula from Jac13
         rho         = 0.29d0
         gam         = 0.689d0
         mus_ref_500 = 66.7d0 
         mus_r       = mus_ref_500 * (wave / 500.d0)**(-4.d0)
         mus_m       = mus_ref_500 * (wave / 500.d0)**(-gam)
         mus         = rho * mus_r + (1.d0 - rho) * mus_m       !in cm-1
         mus         = mus / (1.0d0 - hgg)
         
         Stratum = mua !+ mus
      else
         print*,'Wavelength out of range!'
         print*,'Setting Stratum mua = 0!!!'
         Stratum = 0.d0
      end if
   end function stratum
   
   DOUBLE PRECISION function Epidermis(wave)

      DOUBLE PRECISION, intent(in) :: wave
      DOUBLE PRECISION             :: mua, mus, rho, gam, mus_ref_500, mus_r, mus_m
      DOUBLE PRECISION             :: a, f_ray, bcoeff, nu_m , W, S, B, baseline 
      DOUBLE PRECISION             :: C_bili, C_caro
      

      if(wave .ge. 220.d0)then
      !set mua
      nu_m     = 0.01d0
      W        = 2.0d0      !water level
      S        = 0.75d0     !Blood oxygenation
      B        = 0.0d0      !blood fraction
      C_caro   = 2.1d-4
      C_bili   = 0.0d0
      Baseline = 0.3d0 
                                                                                    
      mua = B * S * Oxy_Hb( wave) + B * (1.d0 - S) * Deoxy_Hb(wave) + W * water(wave) + &
            (7.0d0 - W) * fat(wave) + 2.3d0 * nu_m * (Eumel(wave) + Pheomel(wave)) + &
            2.3d0 * C_bili * bilirubin(wave) + 2.3d0 * C_caro * carotene(wave) + &
            (1.d0 - Baseline) * base(wave)   !in cm-1   
       
      !set mus
      ! formula from Jac13
         rho         = 0.29d0
         gam         = 0.689d0
         mus_ref_500 = 66.7d0 
         mus_r       = mus_ref_500 * (wave / 500.d0)**(-4.d0)
         mus_m       = mus_ref_500 * (wave / 500.d0)**(-gam)
         mus         = rho * mus_r + (1.d0 - rho) * mus_m       !in cm-1
         mus         = mus / (1.0d0 - hgg)
         
         Epidermis = mua! + mus
      else
         print*,'Wavelength out of range!'
         print*,'Setting Epidermis mua = 0!!!'
         Epidermis = 0.d0
      end if
   end function Epidermis
   
   DOUBLE PRECISION function Pap_dermis(wave)

      DOUBLE PRECISION, intent(in) :: wave
      DOUBLE PRECISION             :: mua, mus, a, f_ray, bcoeff,  W, S, B, baseline
      DOUBLE PRECISION             :: C_bili, C_caro, nu_m

      if(wave .ge. 250.d0)then
      !set mua
      nu_m     = 0.0d0
      W        = 5.0d0      !water level
      S        = 0.75d0     !Blood oxygenation
      B        = 6.0d0     !blood fraction
      C_caro   = 7.0d-5
      C_bili   = 0.05d0
      Baseline = 0.3d0
      
      mua = B * S * Oxy_Hb( wave) + B * (1.d0 - S) * Deoxy_Hb(wave) + W * water(wave) + &
            (7.0d0 - W) * fat(wave) + 2.3d0 * nu_m * (Eumel(wave) + Pheomel(wave)) + &
            2.3d0 * C_bili * bilirubin(wave) + 2.3d0 * C_caro * carotene(wave) + &
            (1.d0 - Baseline) * base(wave)   !in cm-1
      !set mus
      ! formula from Jac13                                                                 
         a     = 43.6d0
         f_ray = 0.41d0
         b     = 0.562
         mus   = a * (f_ray * (wave / 500.d0)**(-4.d0) + (1.d0 - f_ray) * (wave / 500.d0)**(-bcoeff))   !in cm-1
         mus   = mus / (1.0d0 - hgg)

         Pap_dermis = mua! + mus
      else
         print*,'Wavelength out of range!'
         print*,'Setting Pap_dermis mua = 0!!!'
         Pap_dermis = 0.d0
      end if
   end function Pap_dermis
   
   DOUBLE PRECISION function Ret_dermis(wave)

      DOUBLE PRECISION, intent(in) :: wave
      DOUBLE PRECISION             :: mua, mus, a, f_ray, bcoeff,  W, S, B, baseline
      DOUBLE PRECISION             :: C_bili, C_caro, nu_m   
      
      if(wave .ge. 250.d0)then
      !set mua
      nu_m     = 0.0d0
      W        = 7.0d0      !water level
      S        = 0.75d0     !Blood oxygenation
      B        = 4.5d0     !blood fraction
      C_caro   = 7.0d-5
      C_bili   = 0.05d0
      Baseline = 0.3d0
      
      mua = B * S * Oxy_Hb( wave) + B * (1.d0 - S) * Deoxy_Hb(wave) + W * water(wave) + &
            (7.0d0 - W) * fat(wave) + 2.3d0 * nu_m * (Eumel(wave) + Pheomel(wave)) + &
            2.3d0 * C_bili * bilirubin(wave) + 2.3d0 * C_caro * carotene(wave) + &
            (1.d0 - Baseline) * base(wave)   !in cm-1
      !set mus
      ! formula from Jac13                                                                 
         a     = 43.6d0
         f_ray = 0.41d0
         b     = 0.562
         mus   = a * (f_ray * (wave / 500.d0)**(-4.d0) + (1.d0 - f_ray) * (wave / 500.d0)**(-bcoeff))         !in cm-1
         mus   = mus / (1.0d0 - hgg)
            
         Ret_dermis = mua!+ mus
      else
         print*,'Wavelength out of range!'
         print*,'Setting Ret_dermis mua = 0!!!'
         Ret_dermis = 0.d0
      end if
   end function Ret_dermis
   
   DOUBLE PRECISION function Hypo_dermis(wave)   
      
      DOUBLE PRECISION, intent(in) :: wave
      DOUBLE PRECISION             :: mua, mus, W, S, B, baseline
      DOUBLE PRECISION             :: C_bili, C_caro, nu_m     
      
      if(wave .ge. 250.d0)then
      !set mua
      nu_m     = 0.0d0
      W        = 7.d0      !water level
      S        = 0.75d0     !Blood oxygenation
      B        = 5.0d0     !blood fraction
      C_caro   = 0.0d0
      C_bili   = 0.0d0
      Baseline = 0.3d0
      
      mua = B * S * Oxy_Hb( wave) + B * (1.d0 - S) * Deoxy_Hb(wave) + W * water(wave) + &
            (7.d0 - W) * fat(wave) + 2.3d0 * nu_m * (Eumel(wave) + Pheomel(wave)) + &
            2.3d0 * C_bili * bilirubin(wave) + 2.3d0 * C_caro * carotene(wave) + &
            (1.d0 - Baseline) * base(wave)   !in cm-1       
      !set mus
         mus = 1050.6d0 * wave**(-0.68d0) !in cm-1
         mus = mus / (1.0d0 - hgg)
         
         Hypo_dermis = mua! + mus
      else
         print*,'Wavelength out of range!'
         print*,'Setting Hypo_dermis mua = 0!!!'
         Hypo_dermis = 0.0d0
      end if
   end function hypo_dermis
end module skin_prop
