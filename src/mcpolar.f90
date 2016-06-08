program mcpolar

!external libs
use mpi

!shared data
use constants
use photon_vars
use iarray
use opt_prop

!subroutines
use subs
use ppm
use reader_mod
use gridset_mod
use sourceph_mod
use noisey_mod
use tauint
use ch_opt
use stokes_mod
use binning_mod
use peelingoff_mod
use writer_mod

implicit none

integer :: nphotons,iseed,j,xcell,ycell,zcell,flucount,nlow
logical :: tflag,sflag,fflag
logical :: fad_bool, nad_bool, nadh_bool, ribo_bool, try_bool, tyro_bool
DOUBLE PRECISION :: ddy,ddx,nscatt,n1,n2,weight,val,fluro_prob
DOUBLE PRECISION :: delta,xcur,ycur,zcur,thetaim,ran, wave_in

DOUBLE PRECISION :: ddz,ddr,v(3)
DOUBLE PRECISION :: costim,cospim,sintim,sinpim
DOUBLE PRECISION :: phiim
real :: start,finish,ran2,sleft,fleft,time

!      variables for openmpi. GLOBAL indicates final values after mpi reduce
!      numproc is number of process being run, id is the indivdual id of each process
!      error is the error flag for MPI

DOUBLE PRECISION :: nscattGLOBAL
integer          :: error,numproc,id

!set directory paths
call directory(id)

call alloc_array(0)
call zarray

!      !init mpi
call MPI_init(error)
!      ! get number of processes
call MPI_Comm_size(MPI_COMM_WORLD,numproc,error)
!      ! get individual process id
call MPI_Comm_rank(MPI_COMM_WORLD,id,error)


!**** Read in parameters from the file input.params
open(10,file=trim(resdir)//'input.params',status='old')
    read(10,*) nphotons
    read(10,*) xmax
    read(10,*) ymax
    read(10,*) zmax
    read(10,*) n1
    read(10,*) n2
    close(10)

! set seed for rnd generator. id to change seed for each process
iseed=95648324+id

!read in optical property data
call reader1
acount=0
fcount=0

!****** setup up arrays and bin numbers/dimensions

!set bin widths for deposit method
ddz=(2.*zmax)/cbinsnum
ddx=(2.*xmax)/cbinsnum
ddy=(2.*ymax)/cbinsnum
ddr=(ymax)/(cbinsnum)

!***** Set up constants, pi and 2*pi  ********************************

iseed=-abs(iseed)  ! Random number seed must be negative for ran2
flucount=0

! postion image in degrees
phiim=0.*pi/180.
thetaim=0.*pi/180.

!     image postion vector
!angle for vector
sintim=sin(thetaim)
sinpim=sin(phiim)
costim=cos(thetaim)
cospim=cos(phiim)

!vector
v(1)=sintim*cospim     
v(2)=sintim*sinpim
v(3)=costim  

!set optical properties and make cdfs.
wave_in = 365.d0 
wave    = wave_in
call opt_set

! calculate number of photons to be run over all cores.  
if(id.eq.0)then
   print*, ''      
   print*,'# of photons to run',nphotons*numproc
end if

!***** Set up density grid *******************************************
call gridset(id, wave)
!call opt_set()
!***** Set small distance for use in optical depth integration routines 
!***** for roundoff effects when crossing cell walls
delta=1.e-6*(2.*xmax/nxg)
nscatt=0
tcount=0

call MPI_Barrier(MPI_COMM_WORLD,error)
call cpu_time(start)
print*, ' '
print*, 'Photons now running on core:',id
call cpu_time(sleft)
call mk_cdf(nadh_fluro,nadh_cdf,size(nadh_fluro,1))
call mk_cdf(fad_fluro,fad_cdf,size(fad_fluro,1)) 
call mk_cdf(ribo_fluro,ribo_cdf,size(ribo_fluro,1)) 
call mk_cdf(try_fluro,try_cdf,size(try_fluro,1)) 
call mk_cdf(tyro_fluro,tyro_cdf,size(tyro_fluro,1)) 


!loop over photons   
do j=1,nphotons
  
!set init weight and flags
   wave=wave_in
   fad_bool=.FALSE.; nad_bool=.FALSE.; nadh_bool=.FALSE.; ribo_bool=.FALSE.; 
   try_bool=.FALSE.; tyro_bool=.FALSE.
   call opt_set
   tflag=.FALSE.
   fflag=.FALSE.
   sflag=.TRUE.     ! flag for fresnel subroutine. so that incoming photons
                       ! are treated diffrently to outgoing ones

! code to output progress as program runs also gives an estimate of when program will complete.
   if(mod(j,int(0.02d0*nphotons)).eq.0)then
      if(id.eq.0)then
          write(*,FMT="(A1,A,t21,F6.2,A,$)") achar(13), &
                " Percent Complete: ", (real(j)/real(nphotons))*100.0, "%"
      end if
   end if
   if(id.eq.0)then
   if (j.eq.1000)then
      call cpu_time(fleft)
      time = ((fleft-sleft)/1000.)*real(nphotons)
      print*,' '
      if(time.ge.60.)then
         print'(A, I3, 1X, A)','Approx time program will take to run: ',floor((time)/60.d0),'mins'
      else
         print'(A, 1X, I2, A)', 'Approx time program will take to run:',floor(time),'s'
      end if
         print*,' '
   end if
   end if
!***** Release photon from point source *******************************
   call sourceph(xcell,ycell,zcell,iseed)

!***** Update xcur etc.
   xcur=xp+xmax
   ycur=yp+ymax
   zcur=zp+zmax

!***** Generate new normal corresponding to bumpy surface
!   call noisey(xcell,ycell)

!***** check whether the photon enters medium     
   call fresnel(n1,n2,sflag,tflag,iseed,ddx,ddy,weight,xcur,ycur)
   sflag=.FALSE.

!****** Find scattering location
   call tauint2(n1,n2,xcell,ycell,zcell,&
   tflag,iseed,delta,sflag,weight,ddx,ddy)
     
!************ Peel off photon into image
!   call peelingoff(xcell,ycell,zcell,delta, &
!   v,sintim,costim,sinpim,cospim)
!   print*,tflag,xp,yp,zp
!******** Photon scatters in grid until it exits (tflag=TRUE) 
   do while(tflag.eqv..FALSE.) 
!******** Scatter or absorb/fluro

! Select albedo based on current photon wavelength
      ran=ran2(iseed)
!      print*,zcell,albedo(xcell,ycell,zcell,1)
      if(ran.lt.albedo(1,1,zcell,1))then !photons scatters
         call stokes(iseed)
         nscatt=nscatt+1
      else !photon absorbs
         if( ran .lt. mu_nadh/rhokap(1, 1, zcell , 1) + albedo(1, 1, zcell, 1))then
!           if(wave.lt.100)print*,'nadh1',wave
            call sample(nadh_fluro,nadh_cdf,wave,iseed)
            call opt_set()
            nadh_bool=.TRUE.
            fad_bool=.FALSE.; nad_bool=.FALSE.; ribo_bool=.FALSE.; 
            try_bool=.FALSE.; tyro_bool=.FALSE.
!            if(wave.lt.100)print*,'nadh2',wave
         elseif( ran .lt. (mu_ribo+mu_nadh)/rhokap(1, 1, zcell , 1) &
                                          + albedo(1, 1, zcell, 1))then
!            if(wave.lt.100)print*,'ribo1',wave
            call sample(ribo_fluro,ribo_cdf,wave,iseed)
            call opt_set()
            ribo_bool=.TRUE.
            fad_bool=.FALSE.; nad_bool=.FALSE.; nadh_bool=.FALSE.; 
            try_bool=.FALSE.; tyro_bool=.FALSE.
!            if(wave.lt.100)print*,'ribo2',wave
          elseif( ran .lt. (mu_ribo+mu_nadh+mu_fad)/rhokap(1, 1, zcell , 1) &
                                          + albedo(1, 1, zcell, 1))then
!            print*,'fad1',wave
            call sample(fad_fluro,fad_cdf,wave,iseed)
            call opt_set()
            fad_bool=.TRUE.
            ribo_bool=.FALSE.; nad_bool=.FALSE.; nadh_bool=.FALSE.; 
            try_bool=.FALSE.; tyro_bool=.FALSE.
!            if(wave.lt.100)print*,'fad2',wave
         elseif( ran .lt. (mu_try+mu_ribo+mu_nadh+mu_fad)/rhokap(1, 1, zcell , 1) &
                                                + albedo(1, 1, zcell, 1))then
!            if(wave.lt.100)print*,'try1',wave
            call sample(try_fluro,try_cdf,wave,iseed)
            call opt_set()
            try_bool=.TRUE.
            fad_bool=.FALSE.; nad_bool=.FALSE.; nadh_bool=.FALSE.; ribo_bool=.FALSE.; 
            tyro_bool=.FALSE.
!           if(wave.lt.100) print*,'try2',wave
         elseif( ran .lt. (mu_tyro+mu_try+mu_ribo+mu_nadh+mu_fad)/rhokap(1, 1, zcell , 1) &
                                                + albedo(1, 1, zcell, 1))then
!            if(wave.lt.100)print*,'tyro',wave
            call sample(tyro_fluro,tyro_cdf,wave,iseed)
            call opt_set()
            tyro_bool=.TRUE.
            fad_bool=.FALSE.; nad_bool=.FALSE.; nadh_bool=.FALSE.; ribo_bool=.FALSE.; 
            try_bool=.FALSE.
!            if(wave.lt.100)print*,'tyro2',wave
         else
            tflag=.TRUE.
         end if
      end if

!maxval(excite_array,2) gives wavelength col
!minval(maxval(excite_array,2)) gives min in wavelength col

!******** Drop weight in appro bin
!      call binning(ddr,zcur,ddz,absorb)

!************ Find next scattering location
      call tauint2(n1,n2,xcell,ycell,zcell,&
   tflag,iseed,delta,sflag,weight,ddx,ddy)

!************ Peel off photon into image
!      call peelingoff(xcell,ycell,zcell,delta &
!      ,v,sintim,costim,sinpim,cospim)

      xcur=xp+xmax
      ycur=yp+ymax
      zcur=zp+zmax

   end do
   if(zp.ge.zmax*.999)then
   flucount=flucount+1
      call wave_to_RGB(wave, xcell, ycell , rgb, 0.8d0)
   end if
!bin photons leaving top surface if they are collected by fibre
   if(int(wave).ne. wave_in .and.nxp.gt. 0.d0)then
      fluroexit(int(wave),7)=fluroexit(int(wave),7)+1
      if(nadh_bool)then
         fluroexit(int(wave),1)=fluroexit(int(wave),1)+1
!         print*,'nadh'
      elseif(ribo_bool)then
         fluroexit(int(wave),2)=fluroexit(int(wave),2)+1
!         print*,'ribo'
      elseif(try_bool)then
         fluroexit(int(wave),3)=fluroexit(int(wave),3)+1
!         print*,'try'   
      elseif(tyro_bool)then
         fluroexit(int(wave),4)=fluroexit(int(wave),4)+1
!         print*,'tyro'
      elseif(fad_bool)then
         fluroexit(int(wave),5)=fluroexit(int(wave),5)+1
!         print*,'fad'
      end if
      flucount=flucount+1
   end if
end do      ! end loop over nph photons
call MPI_BARRIER(MPI_COMM_WORLD, error)
print*,' '
print*,sum(fluroexit(:,7))

print*, flucount,'# photons collected',id

!give time taken to run program
call MPI_BARRIER(MPI_COMM_WORLD, error)
call cpu_time(finish)
if(finish-start.ge.60.)then
 print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
else
      print*, 'time taken ~',floor(finish-start/60.),'s'
end if

!force syncro
call MPI_Barrier(MPI_COMM_WORLD,error)
call alloc_array(1)
rgb(:,:,1)=rgb(:,:,1)/maxval(rgb(:,:,1))
rgb(:,:,2)=rgb(:,:,2)/maxval(rgb(:,:,2))
rgb(:,:,3)=rgb(:,:,3)/maxval(rgb(:,:,3))

!if(id.eq.0)then
!call write_ppm(trim(fileplace)//'render.ppm', rgbGLOBAL)
!end if

!      path length reduce     
call MPI_REDUCE(jmean,jmeanGLOBAL,((nxg+3)*(nyg+3)*(nzg+3))*4,MPI_DOUBLE_PRECISION &
               ,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_Barrier(MPI_COMM_WORLD,error)
print*,'done jmean'

rgb=rgb*256.d0
!call MPI_Allreduce(rgb, rgbGLOBAL, size(rgbGLOBAL), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,error)
print*,'done'

                  
!call MPI_REDUCE(rgb,rgbGLOBAL,size(rgbGLOBAL),MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!               MPI_COMM_WORLD,error)






!call MPI_REDUCE(dep,depGLOBAL,cbinsnum,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
!call MPI_Barrier(MPI_COMM_WORLD,error)

!     deposit reduce
!call MPI_REDUCE(deposit,depositGLOBAL,(cbinsnum**3),MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
!call MPI_Barrier(MPI_COMM_WORLD,error)
!print*,'done deposit'

!!      images reduce
!call MPI_REDUCE(image,imageGLOBAL,(Nbins**2)*4,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
!call MPI_Barrier(MPI_COMM_WORLD,error)
print*,'done image'

!     nscatt reduce
call MPI_Barrier(MPI_COMM_WORLD,error)
call MPI_REDUCE(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)

!     trans reduce
!call MPI_Barrier(MPI_COMM_WORLD,error)
!call MPI_REDUCE(trans,transGLOBAL,nxg*nyg,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
print*,'done trans'

!call MPI_Barrier(MPI_COMM_WORLD,error)
!call MPI_REDUCE(follow,followGLOBAL,(nxg)*(nyg)*(nzg),MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)

!     fluroexit reduce
call MPI_Barrier(MPI_COMM_WORLD,error)
call MPI_REDUCE(fluroexit,fluroexitGLOBAL,1000*7,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)

!     fluro_pos reduce
!call MPI_Barrier(MPI_COMM_WORLD,error)
!call MPI_REDUCE(fluro_pos,fluro_posGLOBAL,nxg*nyg*nzg,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)

print*,'done all reduces',id

call MPI_Barrier(MPI_COMM_WORLD,error)
if(id.eq.0)then
    print*,'Average # of scatters per photon:',sngl(nscattGLOBAL/(nphotons*numproc))

    !write out files

    call writer(nphotons, numproc)
    print*,'write done'
end if

!end MPI processes
call MPI_Finalize(error)
end program mcpolar
