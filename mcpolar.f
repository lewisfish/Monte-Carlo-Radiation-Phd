      program mcpolar

      implicit none

      include 'grid.txt'
      include 'photon.txt'

c***** Parameter declarations ****************************************
      integer nphotons,iseed,j,xcell,ycell,zcell,tflag
      integer cnt,io,sflag,Nbins,cur,cbinsnum,bflag,i
      real*8 nscatt
      real mus1,mua1,hgg,xmax,ymax,zmax,kappa1,albedo,absorb,albedo1
      real pi,twopi,fourpi,g2,delta,xcur,ycur,zcur,d,kappa2
      real, allocatable :: noise(:,:),reflc(:,:),trans(:,:),image(:,:)
      real, allocatable :: flu(:,:),deposit(:,:,:),fluro(:,:,:)
      real ran2,n1,n2,start,finish,weight,terminate,chance,mus2,mua2
      real albedo2,mua,kappa,ddz,ddx,ddy,angle
      
      character(*),parameter::
     +      fileplace="/home/lewis/phdshizz/grid/data/"



! todo 
!add weight,russian roulette, mus/mua in place of kappa/albedo-done
!fix fresnel somehow + bump map
!add kennys stuff(piece of paper) 
! peeling off !!!!!!!!!testing!!!!!!!!!!!!!!!!!!
!add colours/render shizz
!add bins and proper formatting
!convert to f95
!parallize-buy books


c**** Read in parameters from the file input.params
      open(10,file='input.params',status='old')
          read(10,*) nphotons
          read(10,*) iseed
          read(10,*) mua1
          read(10,*) mus1
          read(10,*) hgg
          read(10,*) xmax
          read(10,*) ymax
          read(10,*) zmax
          read(10,*) n1
          read(10,*) n2
          close(10)   


c***** read in noise data

      open(13,file=fileplace//'noisedots.dat')
      cnt=0
      do

        read(13,*,IOSTAT=io)
        
      if (io < 0) then
           close(13)
           allocate(noise(1:cnt,1:cnt))
           noise=0.
           exit
      else
       cnt=cnt+1
      end if

      end do
      open(14,file='noisedots.dat') 
      do i=1,cnt
            read(14,*) (noise(i,j),j=1,cnt)
      end do



c****** setup up arrays and bin numbers/dimensions

      Nbins=401
      cbinsnum=200
     
      ddz=(2.*zmax)/cbinsnum
      ddx=(2.*xmax)/cbinsnum
      ddy=(2.*ymax)/cbinsnum

      allocate(reflc(1:cnt,1:cnt),trans(1:cnt,1:cnt))
      allocate(image(-((Nbins-1)/2):((Nbins-1)/2),-((
     +      Nbins-1)/2):((Nbins-1)/2)))
      allocate(flu(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2)
     +      :((Nbins-1)/2)))
      allocate(deposit(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum))
      allocate(fluro(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum))

      reflc=0.
      trans=0.
      image=0.
      flu=0.
      deposit=0.
      fluro=0.


c***** Set up constants, pi and 2*pi  ********************************
      pi=4.*atan(1.)
      twopi=2.*pi
      fourpi=4.*pi
      sflag=1
      iseed=-abs(iseed)  ! Random number seed must be negative for ran2

      g2=hgg*hgg  ! Henyey-Greenstein parameter, hgg^2
      
      ! set the terminal value for russian roulette and init weight
      terminate=0.0001
      chance=0.1
      mus2=5.
      mua2=.2
      kappa1=mua1+mus1
      kappa2=mua2+mus2
      albedo1=mus1/kappa1
      albedo2=mus2/kappa2
      
c**** Initialize arrays to zero *************************************
      call iarray(xface,yface,zface,rhokap)
      
c***** Set up density grid *******************************************
      call gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,
     +                  kappa1,kappa2)

c***** Set small distance for use in optical depth integration routines 
c***** for roundoff effects when crossing cell walls
      delta=1.e-6*(2.*xmax/nxg)


c**** Loop over nph photons from each source *************************
        nscatt=0
        call cpu_time(start)
        do j=1,nphotons
            !set init weight and optical properties
            weight=1.0
            cur=1
          if(mod(j,10000).eq.0)then
             print *, j,' scattered photons completed'
          end if

c***** Release photon from point source *******************************
          call sourceph(xp,yp,zp,nxp,nyp,nzp,
     +                    sint,cost,sinp,cosp,phi,
     +                    xmax,ymax,zmax,twopi,
     +                    xcell,ycell,zcell,nxg,nyg,nzg,iseed)
     
c***** Update xcur etc.

      xcur=xp+xmax
      ycur=yp+ymax
      zcur=zp+zmax     
     
c***** Generate new normal corresponding to bumpy surface
!          call noisey(xcell,ycell,noise,cnt,angle,nxp,
!     +      nyp,nzp,cost,sint,cosp,sinp)
     
C***** check whether the photon enters medium     
!          call fresnel(n1,n2,cost,sint,sflag,tflag,iseed,
!     +      reflc,xcell,ycell,cnt,trans,weight)

!            sflag=0

c****** Find scattering location
          call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +               kappa2,xface,yface,zface,rhokap,noise,cnt,d,
     +                  xcell,ycell,zcell,tflag,iseed,delta,cur)



c******** Photon scatters in grid until it exits (tflag=1) 
          do while(tflag.eq.0)
          if(rhokap(xcell,ycell,zcell,cur).eq.kappa2)then
                  cur=2
          end if
c******** Drop weight
            if(cur.eq.1)then
                  kappa=kappa1
                  albedo=albedo1
                  mua=mua1
            else
                  kappa=kappa2
                  albedo=albedo2
                  mua=mua2
            end if
            absorb=weight*(mua/kappa)          
            weight=weight*albedo

            bflag=1

            call binning(deposit,bflag,xcur,ycur,weight,ddx,fluro,
     +                  ddy,cbinsnum,zcur,xmax,ymax,zmax,ddz,cur)
          
c************ Scatter photon into new direction and update Stokes parameters
                   call stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
     +                  hgg,g2,pi,twopi,iseed)
     
c************ carry out russian roulette to kill off phototns
            if(weight.le.terminate)then
                  if(ran2(iseed).le.chance)then
                        weight=weight/chance
                  else
                        absorb=weight
           call binning(deposit,bflag,xcur,ycur,weight,ddx,fluro,
     +                  ddy,cbinsnum,zcur,xmax,ymax,zmax,ddz,cur)
                        weight=0.
                        goto 100
                  end if           
            end if                   
                   
            nscatt=nscatt+1

c************ Find next scattering location
              call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +               kappa2,xface,yface,zface,rhokap,noise,cnt,d,
     +                  xcell,ycell,zcell,tflag,iseed,delta,cur)
     
          if(rhokap(xcell,ycell,zcell,cur).eq.kappa2)then
                  cur=2
          end if

          
                   call peelingoff(xmax,ymax,zmax,cur,
     +                  xface,yface,zface,rhokap,
     +                  xcell,ycell,zcell,delta,xp,yp,zp
     +                  ,pi,image,Nbins,flu)
         
            xcur=xp+xmax
            ycur=yp+ymax
            zcur=zp+zmax
            
          end do

100      continue

        end do      ! end loop over nph photons
      call cpu_time(finish)
                  if(finish-start.ge.60.)then
            print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
            else
            print*, 'time taken ~',floor(finish-start/60.),'s'
            end if
      print*,'Avereage number of scatterings = ',sngl(nscatt/nphotons)

      call writer(image,reflc,trans,Nbins,fileplace,flu,deposit,
     +                  cbinsnum,fluro,cnt)

      stop
      end
