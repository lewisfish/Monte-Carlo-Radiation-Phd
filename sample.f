      subroutine opt_sample(fileplace,iseed,output,inputspec,choice)
      
      implicit none
      
      double precision :: value,choice
      double precision, intent(out) :: output
      double precision,allocatable :: wave(:),flux(:),cdf(:)
      double precision :: summ
      real ran2
      integer :: i,j,Nlower
      integer,intent(in) :: iseed
      character(*),intent(in):: fileplace,inputspec

      cdf=0.
      flux=0.
      wave=0.
      summ=0.
      
      open(1,file=fileplace//inputspec)
      cnt=0
      do

        read(1,*,IOSTAT=io)
        
            if (io < 0) then
                  close(1)
                  allocate(wave(cnt),flux(cnt),cdf(cnt-1))
                  exit
            else
                  cnt=cnt+1
            end if

      end do

      open(2,file=fileplace//inputspec)
      do i=1,cnt
            read(2,*) wave(i),flux(i)
      end do

      do i=1,cnt-1
            summ=0.
            do j=1,i
            summ=summ+.5*(flux(j+1)+flux(j))*(wave(j+1)-wave(j))
            end do
            cdf(i)=summ
      end do
      
      cdf=cdf/cdf(cnt-1)
      
      call search(value,cdf,cnt,iseed,Nlower)
      
      if(choice.eq.1.)then
      !select wavelength to fluro at
      value=ran2(iseed)
            call search(value,cdf,cnt,Nlower)
      output=wave(Nlower+1)+(wave(Nlower+2)-wave(Nlower+1))*((value-
     +cdf(Nlower))/(cdf(Nlower+1)-cdf(Nlower)))
      else  !change value of mua
            call search(choice,wave,cnt,Nlower)
      output=flux(Nlower+1)+(flux(Nlower+2)-flux(Nlower+1))*((choice-
     +cdf(Nlower))/(cdf(Nlower+1)-cdf(Nlower)))
      end if

      end program opt_sample
