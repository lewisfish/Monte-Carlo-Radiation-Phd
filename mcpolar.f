      program mcpolar

      use mpi
      implicit none

      include 'grid.txt'
      include 'photon.txt'

c***** Parameter declarations ****************************************
      integer nphotons,iseed,j,xcell,ycell,zcell,ibank
      integer cnt,io,Nbins,cur,cbinsnum,i,flucount
      integer face(6)
      logical tflag,forceflag,sflag,tauflag
      double precision  xmax,ymax,zmax,absorb,dens(8),zlow
      double precision  pi,twopi,fourpi,delta,xcur,ycur,zcur,d,thetaim
      double precision ,allocatable :: noise(:,:),reflc(:,:),trans(:,:)
      double precision ,allocatable :: deposit(:,:,:,:),image(:,:,:)
      double precision ,allocatable :: kappa(:),albedo(:),bank(:,:,:)
      double precision  n1,n2,n3,weight,terminate,nscatt
      double precision  ddz,ddx,ddy,tau,v(3)
      double precision  costim,cospim,sintim,sinpim,xlow,xhi,ylow,yhi
      double precision  mua(8),mus(8),hgg(8),g2(8),phiim,zhi,chance
      real start,finish,ran2
      
!      !variables for openmpi. GLOBAL indicates final values after mpi reduce
!      !numproc is number of process being run, id is the indivdual id of each process
!      !error is the error flag for MPI
      
      double precision nscattGLOBAL
      integer error,numproc,id
      
      double precision, allocatable :: imageGLOBAL(:,:,:)
      double precision, allocatable :: depositGLOBAL(:,:,:,:)
      
      character(len=100) :: opt_parmas
      
!      !set directory for data storage
      character(*),parameter::
     +      fileplace="/home/st-andrews/lm959/data/"
     
      !get filename for optical properties
      call get_command_argument(1,opt_parmas)
      
!      !init mpi
      call MPI_init(error)
!      ! get number of processes
      call MPI_Comm_size(MPI_COMM_WORLD,numproc,error)
!      ! get individual process id
      call MPI_Comm_rank(MPI_COMM_WORLD,id,error)

!!!!!!!!!!!!!!!!!!! todo !!!!!!!!!!!!!!!!
!fix up fluro shizz
!change opt arrays to allocatble 
!create subroutine for reading files, i.e. noise data.
!somehow add all array setting to iarray !!!!!!!!!!!can do this by using modules!!!!!!!!!!!!!!!!
!add in fresnel refections to all surface (for plastic box-interlipid simulation)
!fix fresnel somehow + bump map
!add kennys stuff(piece of paper) 
!add colours/render shizz
!proper formatting
!convert to f95
!parallize-done but not 100% happy with. change makefile so that mpi is an option.

c**** Read in parameters from the file input.params
      open(10,file='input.params',status='old')
          read(10,*) nphotons
          read(10,*) xmax
          read(10,*) ymax
          read(10,*) zmax
          read(10,*) n1
          read(10,*) n2
          read(10,*) n3
          read(10,*) xlow
          read(10,*) xhi
          read(10,*) ylow
          read(10,*) yhi
          read(10,*) zlow
          read(10,*) zhi
          close(10)   
      
      ! set seed for rnd generator. id to change seed for each process
      iseed=95648324+id

      call reader(hgg,mua,mus,opt_parmas,cur)

      dens(1)=1.113
      dens(2)=4.56
      dens(3)=1.113
      dens(4)=4.56
      dens(5)=1.113
      dens(6)=4.56
      dens(7)=1.113
      dens(8)=4.56
c***** read in noise data
      
!      open(13,file=fileplace//'noisedots.dat')
!      cnt=0
!      do

!        read(13,*,IOSTAT=io)
!        
!      if (io < 0) then
!           close(13)
           allocate(noise(1:1,1:1))
!           noise=0.
!           exit
!      else
!       cnt=cnt+1
!      end if

!      end do
!      open(14,file=fileplace//'noisedots.dat') 
!      do i=1,cnt
!            read(14,*) (noise(i,j),j=1,cnt)
!      end do

c****** setup up arrays and bin numbers/dimensions
 
      Nbins=401
      cbinsnum=200
     
!     ! set bin widths for deposit method
      ddz=(2.*zmax)/cbinsnum
      ddx=(2.*xmax)/cbinsnum
      ddy=(2.*ymax)/cbinsnum
      
!      !allocate arrays for program
      allocate(reflc(1:cnt,1:cnt),trans(1:cnt,1:cnt))
      allocate(image(-((Nbins-1)/2):((Nbins-1)/2),-((
     +      Nbins-1)/2):((Nbins-1)/2),4))
      allocate(deposit(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum,4))
      allocate(imageGLOBAL(-((Nbins-1)/2):((Nbins-1)/2),-((
     +      Nbins-1)/2):((Nbins-1)/2),4))
      allocate(depositGLOBAL(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum,4))
      allocate(kappa(cur),albedo(cur))
      allocate(bank(0:10000,1:7,0:1))
!      ! set arrays to 0.
      image=0.
      deposit=0.
      jmean=0.
      jmeanGLOBAL=0.
      imageGLOBAL=0.
      depositGLOBAL=0.
      bank=0.

c***** Set up constants, pi and 2*pi  ********************************
      pi=4.*atan(1.)
      twopi=2.*pi
      fourpi=4.*pi
      sflag=.FALSE.     ! flag for fresnel subroutine. so that incoming photons are treated diffrently to outgoing ones
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
      
c**** Initialize arrays to zero *************************************
      call iarray(xface,yface,zface,rhokap)
      
c***** Set up density grid *******************************************
      call gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,
     +                  kappa,id,cur,xlow,xhi,ylow,yhi,zlow,zhi,face)

c***** Set small distance for use in optical depth integration routines 
c***** for roundoff effects when crossing cell walls
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

          if(mod(j,10000).eq.0)then
             print *, j,' scattered photons completed on core:',id
          end if
          
c***** Release photon from point source if genum is 1*******************************
          call sourceph(xp,yp,zp,nxp,nyp,nzp,
     +                    sint,cost,sinp,cosp,phi,
     +                    xmax,ymax,zmax,twopi,
     +                    xcell,ycell,zcell,nxg,nyg,nzg,iseed)
c***** Update xcur etc.

            xcur=xp+xmax
            ycur=yp+ymax
            zcur=zp+zmax
         
c**** force photon to interact

!      call force(xmax,ymax,zmax,cur,tau,iseed,
!     +                xface,yface,zface,rhokap,nxp,nyp,nzp
!     +                ,xcell,ycell,zcell,delta,xp,yp,zp)           
      forceflag=.FALSE.

c***** Generate new normal corresponding to bumpy surface
!          call noisey(xcell,ycell,noise,cnt,nxp,
!     +      nyp,nzp,cost,sint,cosp,sinp)
     
C***** check whether the photon enters medium     
!          call fresnel(n1,n2,cost,sint,sflag,tflag,iseed,
!     +      reflc,xcell,ycell,cnt,trans,weight,pi)

!            sflag=.FALSE.

c****** Find scattering location
          call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +    forceflag,flucount,kappa,xface,yface,zface,rhokap,noise,
     +    cnt,d,tau,dens,xcell,ycell,zcell,tflag,iseed,jmean,sint,
     +    cost,sinp,phi,twopi,tauflag,weight,pi,sflag,n1,n2,trans,
     +    reflc,cosp,face)
     
c************ Peel off photon into image
          call peelingoff(xmax,ymax,zmax,cur,nxp,nyp,nzp,xface,
     +         yface,zface,rhokap,xcell,ycell,zcell,delta,xp,yp,
     +         zp,pi,image,Nbins,v,g2,hgg,sintim,costim,sinpim,
     +         cospim,tauflag)g)
!      forceflag=.FALSE.

c******** Photon scatters in grid until it exits (tflag=TRUE) 
          do while(tflag.eqv..FALSE.) 
            
c******** Drop weight

      ! Select albedo based on current photon wavelength

            absorb=weight*(mua(cur)/kappa(cur))          
            weight=weight*albedo(cur)
c******** Drop weight in appro bin
            call binning(deposit,xcur,ycur,weight,ddx,
     +                  ddy,cbinsnum,zcur,ddz,cur)
          
c************ Scatter photon into new direction and update Stokes parameters
            if(tauflag.eqv..FALSE.)then
            call stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
     +                  hgg,g2,pi,twopi,iseed,cur)
            end if
c************ carry out russian roulette to kill off phototns
            if(weight.le.terminate)then
                  if(ran2(iseed).le.chance)then
                        weight=weight/chance
                  else
                        absorb=weight
           call binning(deposit,xcur,ycur,weight,ddx,
     +                  ddy,cbinsnum,zcur,ddz,cur)
                        weight=0.
                        exit
                  end if           
            end if                   
                   
            nscatt=nscatt+1

c************ Find next scattering location
          call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +    forceflag,flucount,kappa,xface,yface,zface,rhokap,noise,
     +    cnt,d,tau,dens,xcell,ycell,zcell,tflag,iseed,jmean,sint,
     +    cost,sinp,phi,twopi,tauflag,weight,pi,sflag,n1,n2,trans,
     +    reflc,cosp,face)
     

c************ Peel off photon into image
          call peelingoff(xmax,ymax,zmax,cur,nxp,nyp,nzp,xface,
     +         yface,zface,rhokap,xcell,ycell,zcell,delta,xp,yp,
     +         zp,pi,image,Nbins,v,g2,hgg,sintim,costim,sinpim,
     +         cospim,tauflag)
         
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
      call MPI_REDUCE(jmean,jmeanGLOBAL,((nxg+3)*(nyg+3)*(nzg+3))*4
     +     ,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
     
!     deposit reduce
      call MPI_REDUCE(deposit,depositGLOBAL,(cbinsnum**3)*4,
     +     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)

!      images reduce
      call MPI_REDUCE(image,imageGLOBAL,(Nbins**2)*4,
     +      MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)

!     nscatt reduce
      call MPI_REDUCE(nscatt,nscattGLOBAL,1,
     +    MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
     
      print*,'done all reduces',id
      call MPI_Barrier(MPI_COMM_WORLD,error)
      if(id.eq.0)then
      print*,'Average # of scatters per photon:',sngl(nscattGLOBAL
     +                                        /(nphotons*numproc))

      !write out files
      
      call writer(imageGLOBAL,reflc,trans,Nbins,fileplace,
     +            depositGLOBAL,cbinsnum,cnt,jmeanGLOBAL)
      print*,'write done'
      end if
      
      !end MPI processes
      call MPI_Finalize(error)
      end program mcpolar
