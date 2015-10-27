      subroutine change_opt_prop(optflag)
      
      implicit none
      
      real,intent(in) :: hgg
      real,intent(out) :: mus,mua
      real :: a,b,blanc
      integer,intent(in) :: iseed
      integer,intent(out) ::
      integer :: 
      character(*),intent(in) :: spec1,fileplace
         
      !is g wavelength dep???
      mus=(a*(lam/500.)**(-b))/(1.-hgg)
      
      call opt_sample(fileplace,iseed,mua,'mua.dat',lam)
      
      call update_prop()
      
      end subroutine
