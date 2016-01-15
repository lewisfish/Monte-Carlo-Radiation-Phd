program mcpolar

use mpi
use constants
use photon
use iarray
use subs

implicit none



!***** Parameter declarations ****************************************
integer nphotons,iseed,j,xcell,ycell,zcell,celli,cellk
integer cnt,io,cur,i,flucount,k
logical tflag,sflag,tauflag
DOUBLE PRECISION nscatt
real  xmax,ymax,zmax,absorb,ddy,ddx
real  delta,xcur,ycur,zcur,d,thetaim
real  n1,n2,weight,terminate

real  ddz,ddr,tau,v(3)
real  costim,cospim,sintim,sinpim
real  mua(1),mus(1),hgg(1),g2(1),phiim,chance
real start,finish,ran2

!      !variables for openmpi. GLOBAL indicates final values after mpi reduce
!      !numproc is number of process being run, id is the indivdual id of each process
!      !error is the error flag for MPI

DOUBLE PRECISION nscattGLOBAL
integer error,numproc,id

!set directory paths
call directory

call alloc_array
call zarray

!      !init mpi
call MPI_init(error)
!      ! get number of processes
call MPI_Comm_size(MPI_COMM_WORLD,numproc,error)
!      ! get individual process id
call MPI_Comm_rank(MPI_COMM_WORLD,id,error)

!!!!!!!!!!!!!!!!!!! todo !!!!!!!!!!!!!!!!
!sort fresnel so can differ between transmission on diffrent faces
!investigate +1 in celli,j,k calc in tauint2
!fix up fluro shizz
!change opt arrays to allocatble 
!create subroutine for reading files, i.e. noise data.
!test fresnel
!test bump map, set cross as noise shape then look at diffrent slices/ do calculations manualy
!add kennys stuff(piece of paper) 
!add colours/render shizz
!proper formatting
!parallize-done but not 100% happy with. change makefile so that mpi is an option.

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

call reader(hgg,mua,mus,cur)

!***** read in noise data

open(13,file=trim(resdir)//'noisedots.dat')
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
open(14,file=trim(resdir)//'noisedots.dat') 
do i=1,cnt
read(14,*) (noise(i,j),j=1,cnt)
end do

!****** setup up arrays and bin numbers/dimensions
     
!     ! set bin widths for deposit method
ddz=(2.*zmax)/cbinsnum
ddx=(2.*xmax)/cbinsnum
ddy=(2.*ymax)/cbinsnum
ddr=(ymax+xmax)/(2.*cbinsnum)


!      allocate arrays for program


allocate(kappa(1),albedo(1))

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

! set the terminal value for russian roulette and init weight
terminate=0.0001
! value for russian roulette to increase packet 'energy' by 1/chance      
chance=0.1

!set optical properties
do i=1,cur
g2(i)=hgg(i)**2
kappa(i)=mua(i)+mus(i)
albedo(i)=mus(i)/kappa(i)

end do

if(id.eq.0)then
print*, ''      
print*,'# of photons to run',nphotons*numproc
do i=1,cur
            write(*,*) mua(i),mus(i)
end do
end if

!***** Set up density grid *******************************************
call gridset(xmax,ymax,zmax,kappa,id,cur)

!***** Set small distance for use in optical depth integration routines 
!***** for roundoff effects when crossing cell walls
delta=1.e-6*(2.*xmax/nxg)
nscatt=0

call MPI_Barrier(MPI_COMM_WORLD,error)
call cpu_time(start)
print*, ' '
print*, 'Photons now running on core:',id
  
!loop over photons   
do j=1,nphotons
  
!set init weight and optical properties of current medium(usually 1)

      weight=1.0
      cur=1
      tauflag=.FALSE.
      tflag=.FALSE.
      sflag=.TRUE.     ! flag for fresnel subroutine. so that incoming photons
                       ! are treated diffrently to outgoing ones

    if(mod(j,10000).eq.0)then
 print *, j,' scattered photons completed on core:',id
    end if
    
!***** Release photon from point source *******************************
    call sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed)
!***** Update xcur etc.

xcur=xp+xmax
ycur=yp+ymax
zcur=zp+zmax

!***** Generate new normal corresponding to bumpy surface
    call noisey(xcell,ycell,cnt)

!***** check whether the photon enters medium     
    call fresnel(n1,n2,sflag,tflag,iseed,ddx,ddy,weight,xcur,ycur)
sflag=.FALSE.

!****** Find scattering location
      call tauint2(xmax,ymax,zmax,n1,n2,xcell,ycell,zcell,&
tflag,iseed,delta,sflag,weight,ddx,ddy,cnt)
     
!************ Peel off photon into image
      call peelingoff(xmax,ymax,zmax,xcell,ycell,zcell,delta, &
      v,g2,hgg,sintim,costim,sinpim,cospim)

!******** Photon scatters in grid until it exits (tflag=TRUE) 
    do while(tflag.eqv..FALSE.) 

!******** Drop weight

! Select albedo based on current photon wavelength

absorb=weight*(mua(1)/kappa(1))          
weight=weight*albedo(1)
!******** Drop weight in appro bin
call binning(ddr,zcur,ddz,absorb)
    
!************ Scatter photon into new direction and update Stokes parameters

call stokes(hgg,g2,iseed,cur)
!************ carry out russian roulette to kill off phototns
if(weight.le.terminate)then
      if(ran2(iseed).le.chance)then
            weight=weight/chance
      else
            absorb=weight
     call binning(ddr,zcur,ddz,absorb)
            weight=0.
            tflag=.TRUE.
            exit
      end if           
end if                   
       
nscatt=nscatt+1

!************ Find next scattering location
    call tauint2(xmax,ymax,zmax,n1,n2,xcell,ycell,zcell &
    ,tflag,iseed,delta,sflag,weight,ddx,ddy,cnt)
     

!************ Peel off photon into image
    call peelingoff(xmax,ymax,zmax,xcell,ycell,zcell,delta &
    ,v,g2,hgg,sintim,costim,sinpim,cospim)
   
xcur=xp+xmax
ycur=yp+ymax
zcur=zp+zmax

end do

continue

end do      ! end loop over nph photons
 
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
do i=1,cbinsnum
      dep(k)=dep(k)+deposit(i,k)
end do
end do

dep=dep/(mua(1)*nphotons*(ddz))

call MPI_REDUCE(dep,depGLOBAL,cbinsnum,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_Barrier(MPI_COMM_WORLD,error)

!     deposit reduce
call MPI_REDUCE(deposit,depositGLOBAL,(cbinsnum**2),MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,error)
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


print*,'done all reduces',id
call MPI_Barrier(MPI_COMM_WORLD,error)
if(id.eq.0)then
print*,'Average # of scatters per photon:',sngl(nscattGLOBAL/(nphotons*numproc))

!write out files

call writer
print*,'write done'
end if

!end MPI processes
call MPI_Finalize(error)
end program mcpolar
