MODULE skin_prop

use ch_opt
use opt_prop, only : hgg

implicit none

CONTAINS
   DOUBLE PRECISION function stratum(wave)
   
      DOUBLE PRECISION, intent(in) :: wave
      DOUBLE PRECISION             :: mua, rho, gam, mus_ref_500, mus_r, mus_m, mus


      if(wave.gt.0.0d0)then
      !set mua
      ! formula from GAJG15 supp. material
         mua = ((0.1 - 0.3*10**(-4.) * wave) + 0.125 * (wave/10.) * Base(wave)) * (1. - 0.05) + water(wave)   
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
      DOUBLE PRECISION             :: nu_m, nu_Pd_Hb, nu_Rd_Hb   
      
      nu_m     = 0.001d0
      if(wave .gt. 220.d0)then
      !set mua                                                                                
         mua = (nu_m * (Eumel(wave) + Pheomel(wave)) + (1.d0 - nu_m) * (Carotene(2.1d-4, wave) &
               + (1.d0 - 2.1d-4) * base(wave))) * (1.d0 - .2) + water(wave)             
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
      DOUBLE PRECISION             :: mua, mus, a, f_ray, b, nu_Pd_Hb

      if(wave .gt. 250.d0)then
      !set mua
      ! formula from GAJG15 supp. material
         nu_Pd_Hb = 6.0d0
         mua = (nu_Pd_Hb * (Oxy_Hb(1.d0, wave) + Deoxy_Hb(1.d0, wave) + Bilirubin(7.d0 * 10**(-5), wave) + &
               Carotene(7.d0 * 10**(-5), wave) + (1.- .15) * base(wave))) * (1. - .5) + water(wave) 
      !set mus
      ! formula from Jac13                                                                 
         a     = 43.6d0
         f_ray = 0.41d0
         b     = 0.562
         mus   = a * (f_ray * (wave / 500.d0)**(-4.d0) + (1.d0 - f_ray) * (wave / 500.d0)**(-b))
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
      DOUBLE PRECISION             :: mua, mus, a, f_ray, b, nu_Rd_hb
      
      
      if(wave .gt. 250.d0)then
      !set mua
      ! formula from GAJG15 supp. material
         nu_Rd_hb = 4.5d0
         mua = (nu_Rd_Hb * (Oxy_Hb(1.d0, wave) + Deoxy_Hb(1.d0, wave) + Bilirubin(7.d0*10**(-5), wave) + &
               Carotene(7.d0*10**(-5), wave) + (1.-.15) * base(wave))) * (1.-.7) + water(wave)
      !set mus
      ! formula from Jac13                                                                 
         a     = 43.6d0
         f_ray = 0.41d0
         b     = 0.562
         mus   = a * (f_ray * (wave / 500.d0)**(-4.d0) + (1.d0 - f_ray) * (wave / 500.d0)**(-b))
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
      DOUBLE PRECISION             :: mua, mus, nu_H2O, nu_b, S, b_frac
      
      
      if(wave .gt. 250.d0)then
      !set mua
      nu_H2O = 0.7d0
      nu_b   = 0.05d0
      S      = 0.75d0
      b_frac = 0.002d0 / nu_b
      
      mua = (1.d0 - S) * b_frac * Deoxy_Hb(1.d0, wave) + S * b_frac * nu_b * Oxy_Hb(1.d0, wave) &
             + (1.d0 - b_frac * nu_b) * nu_H2O * water(wave) + (1.d0 - b_frac * nu_b) * &
               (1.d0 - nu_H2O) * base(wave)
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
