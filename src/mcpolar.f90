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

integer nphotons,iseed,j,xcell,ycell,zcell,celli,cellk
integer cnt,io,i,flucount,k,nlow
logical tflag,sflag,fflag
DOUBLE PRECISION nscatt
real :: absorb,ddy,ddx
real :: delta,xcur,ycur,zcur,thetaim,ran
real :: n1,n2,weight,hggtmp

real :: ddz,ddr,v(3),fluro_prob
real :: costim,cospim,sintim,sinpim
real :: phiim
real :: start,finish,ran2

!      variables for openmpi. GLOBAL indicates final values after mpi reduce
!      numproc is number of process being run, id is the indivdual id of each process
!      error is the error flag for MPI

DOUBLE PRECISION :: nscattGLOBAL
integer          :: error,numproc,id

!set directory paths
call directory(id)

call alloc_array
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

!     set bin widths for deposit method
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
wave=405. ! value used in http://www.photobiology.com/v1/sikorski/ and where the excite and fluro data comes from.

call mk_cdf(fluro_array,f_cdf,size(f_cdf))
call mk_cdf(excite_array,e_cdf,size(e_cdf))
call init_opt
   
if(id.eq.0)then
   print*, ''      
   print*,'# of photons to run',nphotons*numproc
end if

!***** Set up density grid *******************************************
call gridset(id)

!***** Set small distance for use in optical depth integration routines 
!***** for roundoff effects when crossing cell walls
delta=1.e-6*(2.*xmax/nxg)
nscatt=0
tcount=0;tcount=0
call MPI_Barrier(MPI_COMM_WORLD,error)
call cpu_time(start)
print*, ' '
print*, 'Photons now running on core:',id

!loop over photons   
do j=1,nphotons
  
!set init weight and flags
   wave=405.
   call init_opt
   tflag=.FALSE.
   fflag=.FALSE.
   sflag=.TRUE.     ! flag for fresnel subroutine. so that incoming photons
                       ! are treated diffrently to outgoing ones

   if(mod(j,100000).eq.0)then
      print *, j,' scattered photons completed on core:',id
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

!******** Photon scatters in grid until it exits (tflag=TRUE) 
   do while(tflag.eqv..FALSE.) 

!******** Drop weight

! Select albedo based on current photon wavelength
      ran=ran2(iseed)
      if(ran.lt.albedo)then !photons scatters
         call stokes(iseed)
         nscatt=nscatt+1
         acount=acount+1
!         print*,mua,mus,kappa,wave
!      else if(ran.lt.(mua+mus)/kappa)then  !photon fluros
       elseif(fflag.eqv..FALSE..and.ran.lt.(mus+mua)/kappa)then
         fflag=.TRUE.
         call sample(excite_array,size(e_cdf),e_cdf,wave,iseed)
         call init_opt
         fluro_pos(xcell,ycell,zcell)=fluro_pos(xcell,ycell,zcell)+1
         fcount=fcount+1
         cost=2.*ran2(iseed)-1.
         sint=(1.-cost*cost)
         if(sint.le.0.)then
            sint=0.
         else
            sint=sqrt(sint)
         endif

         phi=twopi*ran2(iseed)
         cosp=cos(phi)
         sinp=sin(phi)

         nxp=sint*cosp  
         nyp=sint*sinp
         nzp=cost
      else !photon absorbs
!         acount=acount+1
         tflag=.TRUE.
!         
         
!maxval(excite_array,2) gives wavelength col
!minval(maxval(excite_array,2)) gives min in wavelength col

      end if
!******** Drop weight in appro bin
!      call binning(ddr,zcur,ddz,absorb)


!nscatt=nscatt+1
!      print*,tflag,'bt'
!************ Find next scattering location
      call tauint2(n1,n2,xcell,ycell,zcell,&
   tflag,iseed,delta,sflag,weight,ddx,ddy)
!      print*,tflag,'at'

!************ Peel off photon into image
!      call peelingoff(xcell,ycell,zcell,delta &
!      ,v,sintim,costim,sinpim,cospim)
      xcur=xp+xmax
      ycur=yp+ymax
      zcur=zp+zmax

   end do

   continue
   if(int(wave).ne.405..and.zp.ge.zmax-.1)then
      fluroexit(int(wave))=fluroexit(int(wave))+1
   end if
end do      ! end loop over nph photons
print*, acount, '# absorbed',id
print*, fcount, '# fluro',id
print*, acount+fcount,'total',id
call cpu_time(finish)
if(finish-start.ge.60.)then
 print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
else
      print*, 'time taken ~',floor(finish-start/60.),'s'
end if

!force syncro
call MPI_Barrier(MPI_COMM_WORLD,error)

!      path length reduce     
call MPI_REDUCE(jmean,jmeanGLOBAL,((nxg+3)*(nyg+3)*(nzg+3))*4,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_Barrier(MPI_COMM_WORLD,error)
print*,'done jmean'

do k=1,cbinsnum
    do j=1,cbinsnum
        do i=1,cbinsnum
            dep(k)=dep(k)+deposit(i,j,k)
        end do
    end do
end do

!dep=dep/(mua*nphotons*(ddz))

call MPI_REDUCE(dep,depGLOBAL,cbinsnum,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_Barrier(MPI_COMM_WORLD,error)

!     deposit reduce
call MPI_REDUCE(deposit,depositGLOBAL,(cbinsnum**3),MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_Barrier(MPI_COMM_WORLD,error)
print*,'done deposit'

!      images reduce
call MPI_REDUCE(image,imageGLOBAL,(Nbins**2)*4,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_Barrier(MPI_COMM_WORLD,error)
print*,'done image'

!     nscatt reduce
call MPI_Barrier(MPI_COMM_WORLD,error)
call MPI_REDUCE(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)

!     trans reduce
call MPI_Barrier(MPI_COMM_WORLD,error)
call MPI_REDUCE(trans,transGLOBAL,nxg*nyg,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
print*,'done trans'
!     fluroexit reduce
call MPI_Barrier(MPI_COMM_WORLD,error)
call MPI_REDUCE(fluroexit,fluroexitGLOBAL,1000,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
!     fluro_pos reduce
call MPI_Barrier(MPI_COMM_WORLD,error)
call MPI_REDUCE(fluro_pos,fluro_posGLOBAL,nxg*nyg*nzg,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)

print*,'done all reduces',id
call MPI_Barrier(MPI_COMM_WORLD,error)
if(id.eq.0)then
    print*,'Average # of scatters per photon:',sngl(nscattGLOBAL/(nphotons*numproc))

    !write out files

    call writer
    print*,'write done'
end if
print*,'t',tcount,id
print*,'b',bcount,id
!end MPI processes
call MPI_Finalize(error)
end program mcpolar
