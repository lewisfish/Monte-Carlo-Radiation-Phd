
!!!!!!!!!!!!!!!!!!! todo !!!!!!!!!!!!!!!!
!investigate jmean vs deposit method i.e. find normilisation for jmean under russian roulette.
! change reader1 to be better i.e. no repeated code.
!sort fresnel so can differ between transmission on diffrent faces
!investigate +1 in celli,j,k calc in tauint2
!test fresnel
!test bump map, set cross as noise shape then look at diffrent slices/ do calculations manualy
!add kennys stuff(piece of paper) 
!add colours/render shizz
!proper formatting
!parallize-done but not 100% happy with. change makefile so that mpi is an option. ffp.

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
logical tflag,sflag,isoflag
DOUBLE PRECISION nscatt
real :: xmax,ymax,zmax,absorb,ddy,ddx
real :: delta,xcur,ycur,zcur,thetaim,ran
real :: weight,terminate,hggtmp,fdep,hgghold

real :: ddz,ddr,v(3),fluro_prob
real :: costim,cospim,sintim,sinpim
real :: phiim,chance
real :: start,finish,ran2

!      variables for openmpi. GLOBAL indicates final values after mpi reduce
!      numproc is number of process being run, id is the indivdual id of each process
!      error is the error flag for MPI

DOUBLE PRECISION :: nscattGLOBAL
integer          :: error,numproc,id,counter1,counter2

!set directory paths
call directory

!allocate and set arrays to 0
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
!call reader1

counter1=0
counter2=0


!****** setup up arrays and bin numbers/dimensions

!     ! set bin widths for deposit method
ddz=(2.*zmax)/cbinsnum
ddx=(2.*xmax)/cbinsnum
ddy=(2.*ymax)/cbinsnum
ddr=(ymax)/(cbinsnum)

!***** Set up constants, pi and 2*pi  ********************************

iseed=-abs(iseed)  ! Random number seed must be negative for ran2
flucount=0

call init_opt2
   
if(id.eq.0)then
   print*, ''      
   print*,'# of photons to run',nphotons*numproc
end if

!***** Set up density grid *******************************************
call gridset(xmax,ymax,zmax,id)

!***** Set small distance for use in optical depth integration routines 
!***** for roundoff effects when crossing cell walls
delta=1.e-6*(2.*xmax/nxg)

nscatt=0

call MPI_Barrier(MPI_COMM_WORLD,error)
call cpu_time(start)

print*, ' '
print*, 'Photons now running on core:',id
print*,

!loop over photons 
!do nlow=1,2
! jmean=0. 
do j=1,nphotons
  
!set init weight and flags
   !set photon depth flag to -99 so not to mix it up and count normal photons
   fdep=-99.
!   if (nlow.eq.1)then
!      call init_opt1
!   else
!      call init_opt2
!   end if
   call init_opt2
   tflag=.FALSE.
   sflag=.TRUE.     ! flag for fresnel subroutine. so that incoming photons
                    ! are treated diffrently to outgoing ones

   if(mod(j,100000).eq.0)then
      print *, j,' scattered photons completed on core:',id
   end if
    
!***** Release photon from point source *******************************
   call sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed)
!***** Update xcur etc.

   xcur=xp+xmax
   ycur=yp+ymax
   zcur=zp+zmax


!***** check whether the photon enters medium     
   call fresnel(sflag,tflag,iseed,weight,xcur,ycur,xmax,ymax &
               ,zmax,delta,xcell,ycell,zcell)
   sflag=.FALSE.

!****** Find scattering location
   call tauint2(xmax,ymax,zmax,xcell,ycell,zcell,&
   tflag,iseed,delta,sflag,weight,1)


!******** Photon scatters in grid until it exits (tflag=TRUE) 
   do while(tflag.eqv..FALSE.) 
      ran = ran2(iseed)
      if(ran .lt. albedo)then
         call stokes(iseed)
         nscatt = nscatt + 1
      else
         if(fdep .eq. -99)then
            if(ran.lt.(mus + 5.3d-3) / kappa)then
               counter2 = counter2 + 1
               fluro_pos(xcell,ycell,zcell) = fluro_pos(xcell,ycell,zcell) + 1
               fdep = zcell
               depthall(int(fdep)) = depthall(int(fdep)) + 1
              
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
               
               mus = 21. / (1. - hgg)
               mua = 0.23
               kappa = mua + mus + 5.3d-3
               albedo = mus / kappa
            else
               tflag=.true.
               exit
            end if
         else
            tflag=.TRUE.
            exit
         end if
      end if

!************ Find next scattering location
      call tauint2(xmax,ymax,zmax,xcell,ycell,zcell &
      ,tflag,iseed,delta,sflag,weight,1)
      

   end do
   if(fdep.ne.-99.and.zp.ge.zmax-.06)then
      depth(int(fdep))=depth(int(fdep))+1
   end if
end do      ! end loop over nph photons

!end do
print*, counter2, '# fluro',id,sum(depthall),sum(depth)
print*, counter1+counter2,'total',id
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


call MPI_REDUCE(dep,depGLOBAL,cbinsnum,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_Barrier(MPI_COMM_WORLD,error)

call MPI_REDUCE(depth,depthGLOBAL,nzg,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_Barrier(MPI_COMM_WORLD,error)

call MPI_REDUCE(depthall,depthallGLOBAL,nzg,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
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
    call writer(xmax,ymax,zmax,nphotons,numproc)
    print*,'write done'
end if

!end MPI processes
call MPI_Finalize(error)
end program mcpolar
