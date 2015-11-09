      program mcpolar

      use mpi
      implicit none

      include 'grid.txt'
      include 'photon.txt'

c***** Parameter declarations ****************************************
      integer nphotons,iseed,j,xcell,ycell,zcell,ibank,gennum,deathexit
      integer cnt,io,Nbins,cur,cbinsnum,i,flucount,totcount,totphot
      integer face(6)
      double precision  nscatt
      logical tflag,forceflag,sflag,tauflag,stretchflag,loopflag
      double precision  xmax,ymax,zmax,absorb,dens(8),sfact,p,zlow
      double precision  pi,twopi,fourpi,delta,xcur,ycur,zcur,d,thetaim
      double precision , allocatable :: noise(:,:),reflc(:,:),trans(:,:)
      double precision , allocatable :: flu(:,:,:),deposit(:,:,:)
      double precision , allocatable :: kappa(:),albedo(:),bank(:,:,:)
      double precision , allocatable :: image(:,:),fluro(:,:,:,:)
      double precision  n1,n2,n3,weight,terminate
      double precision  ddz,ddx,ddy,tau,v(3),WH,WL,xexit,yexit,zexit
      double precision  costim,cospim,sintim,sinpim,xlow,xhi,ylow,yhi
      double precision  mua(8),mus(8),hgg(8),g2(8),phiim,zhi,chance
      real start,finish,ran2
      
!      !variables for openmpi. GLOBAL indicates final values after mpi reduce
!      !numproc is number of process being run, id is the indivdual id of each process
!      !error is the error flag for MPI
      
      double precision nscattGLOBAL
      double precision flucountGLOBAL,xexitGLOBAL,yexitGLOBAL
      integer error,numproc,id,deathexitGLOBAL,zexitGLOBAL
      
      double precision, allocatable :: imageGLOBAL(:,:),fluGLOBAL(:,:,:)
      double precision, allocatable :: fluroGLOBAL(:,:,:,:)
      double precision, allocatable :: depositGLOBAL(:,:,:)
      
      character(len=100) :: opt_parmas
      
!      !set directory for data storage
      character(*),parameter::
     +      fileplace="/home/lewis/phdshizz/grid/data/"
     
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
!convert real to real*8(or equv)
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
     +      Nbins-1)/2):((Nbins-1)/2)))
      allocate(flu(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2)
     +      :((Nbins-1)/2),4))
      allocate(deposit(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum))
      allocate(fluro(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum,4))

      allocate(imageGLOBAL(-((Nbins-1)/2):((Nbins-1)/2),-((
     +      Nbins-1)/2):((Nbins-1)/2)))
      allocate(fluGLOBAL(-((Nbins-1)/2):((Nbins-1)/2),-((Nbins-1)/2)
     +      :((Nbins-1)/2),4))
      allocate(depositGLOBAL(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum))
      allocate(fluroGLOBAL(-1:cbinsnum,-1:cbinsnum,-1:cbinsnum,4))
      allocate(kappa(cur),albedo(cur))
      allocate(bank(0:10000,1:7,0:1))
!      ! set arrays to 0.
      image=0.
      flu=0.
      deposit=0.
      fluro=0.
      jmean=0.
      jmeanGLOBAL=0.
      imageGLOBAL=0.
      fluGLOBAL=0.
      depositGLOBAL=0.
      fluroGLOBAL=0.
      flucountGLOBAL=0.
      bank=0.
      xexit=0
      yexit=0
      zexit=0
      xexitGLOBAL=0
      yexitGLOBAL=0
      zexitGLOBAL=0
      deathexitGLOBAL=0
      deathexit=0

c***** Set up constants, pi and 2*pi  ********************************
      pi=4.*atan(1.)
      twopi=2.*pi
      fourpi=4.*pi
      sflag=.FALSE.     ! flag for fresnel subroutine. so that incoming photons are treated diffrently to outgoing ones
      iseed=-abs(iseed)  ! Random number seed must be negative for ran2
      flucount=0
      !weight window parameters
      WH=5.
      WL=0.001
      totcount=1
      totphot=0
      loopflag=.FALSE.
      gennum=1
      !exp trasform parameter
      p=0.0
      
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
        call cpu_time(start)
              print*, ' '
              print*, 'Photons now running on core:',id
              
              
c*********************************************************************
c*********************************************************************
c*********************************************************************              
c**** loop over generations of banked photons
        do i=1,gennum
            
            !change bank from one to another when changing generations
!            if(mod(gennum,2).eq.0)then
!                  ibank=1
!            else
!                  ibank=0
!            end if
!            
!            if(i.gt.1)then
!            print*,'going to run:',totcount,'photons. On core:',id
!            totcount=1
!            end if
        !loop over photons   
        do j=1,nphotons
        
            !set init weight and optical properties of current medium(usually 1)
            if(loopflag.eqv..FALSE.)then
                  weight=1.0
                  cur=1
                  tauflag=.FALSE.
            end if

          if(mod(j,10000).eq.0)then
             print *, j,' scattered photons completed on core:',id
          end if
          
!      if(i.eq.1)then
c***** Release photon from point source if genum is 1*******************************
          call sourceph(xp,yp,zp,nxp,nyp,nzp,
     +                    sint,cost,sinp,cosp,phi,
     +                    xmax,ymax,zmax,twopi,
     +                    xcell,ycell,zcell,nxg,nyg,nzg,iseed)
c***** Update xcur etc.

            xcur=xp+xmax
            ycur=yp+ymax
            zcur=zp+zmax
!      else
!      ! else release from bank
!      ! need to change ibank value so to release banked photons from last gen
!            if(ibank.eq.0)ibank=1
!            if(ibank.eq.1)ibank=0
!            xcur=bank(j,1,ibank)
!            ycur=bank(j,2,ibank)
!            zcur=bank(j,3,ibank)
!            
!            nxp=bank(j,4,ibank)
!            nyp=bank(j,5,ibank)
!            nzp=bank(j,6,ibank)
!            
!            weight=bank(j,7,ibank)
!      
!            xp=xcur-xmax
!            yp=ycur-ymax
!            zp=zcur-zmax
!      
!      end if
      
      
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
          call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,forceflag,
     +    flucount,kappa,xface,yface,zface,rhokap,noise,cnt,d,tau,dens,
     +    p,sfact,xcell,ycell,zcell,tflag,iseed,delta,cur,jmean,
     +    stretchflag,sint,cost,sinp,phi,twopi,tauflag,weight,xexit,
     +    yexit,zexit,pi,sflag,n1,n2,trans,reflc,cosp,face)
     
c************ Peel off photon into image
                   call peelingoff(xmax,ymax,zmax,cur,nxp,nyp,nzp
     +                  ,xface,yface,zface,rhokap,
     +                  xcell,ycell,zcell,delta,xp,yp,zp,
     +                  pi,image,Nbins,flu,v,g2,hgg,
     +                  sintim,costim,sinpim,cospim,tauflag)
!      forceflag=.FALSE.


c******** Photon scatters in grid until it exits (tflag=TRUE) 
          do while(tflag.eqv..FALSE.) 
            
c******** Drop weight

      ! Select albedo based on current photon wavelength
            if(stretchflag.eqv..TRUE.)then

            weight=weight*((exp(-kappa(cur)*d*sfact))/(1.-sfact))
            weight=weight*albedo(cur)
            stretchflag=.FALSE.
            !window splitting routine
!            call banking(bank,xcur,ycur,zcur,weight,nxp,nyp,nzp
!     +                  ,WH,WL,totcount,ibank)
            else
            absorb=weight*(mua(cur)/kappa(cur))          
            weight=weight*albedo(cur)
            end if
c******** Drop weight in appro bin
            call binning(deposit,xcur,ycur,weight,ddx,fluro,
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
           call binning(deposit,xcur,ycur,weight,ddx,fluro,
     +                  ddy,cbinsnum,zcur,ddz,cur)
                        weight=0.
                        deathexit=deathexit+1
                        exit
                  end if           
            end if                   
                   
            nscatt=nscatt+1

c************ Find next scattering location
          call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,forceflag,
     +    flucount,kappa,xface,yface,zface,rhokap,noise,cnt,d,tau,dens,
     +    p,sfact,xcell,ycell,zcell,tflag,iseed,delta,cur,jmean,
     +    stretchflag,sint,cost,sinp,phi,twopi,tauflag,weight,xexit,
     +    yexit,zexit,pi,sflag,n1,n2,trans,reflc,cosp,face)
     

c************ Peel off photon into image
                   call peelingoff(xmax,ymax,zmax,cur,nxp,nyp,nzp
     +                  ,xface,yface,zface,rhokap,
     +                  xcell,ycell,zcell,delta,xp,yp,zp,
     +                  pi,image,Nbins,flu,v,g2,hgg,
     +                  sintim,costim,sinpim,cospim,tauflag)
         
            xcur=xp+xmax
            ycur=yp+ymax
            zcur=zp+zmax
            
          end do

      continue

        end do      ! end loop over nph photons
              totphot=totphot+nphotons
              nphotons=totcount
              loopflag=.TRUE.
             
        end do      ! end loop over generations
        
        
      call cpu_time(finish)
            if(finish-start.ge.60.)then
             print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
            else
                  print*, 'time taken ~',floor(finish-start/60.),'s'
            end if


!      path length reduce     
      call MPI_REDUCE(jmean,jmeanGLOBAL,((nxg+3)*(nyg+3)*(nzg+3))*4
     +     ,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)

      call MPI_REDUCE(deposit,depositGLOBAL,cbinsnum**3,
     +     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
      call MPI_REDUCE(fluro,fluroGLOBAL,(cbinsnum**3)*4,
     +     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
     
!      images reduce
      call MPI_REDUCE(image,imageGLOBAL,Nbins**2,
     +      MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
      call MPI_REDUCE(flu,fluGLOBAL,(Nbins**2)*4,
     +     MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
     
!     nscatt reduce
      call MPI_REDUCE(nscatt,nscattGLOBAL,1,
     +    MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
     
!     exit location counter reduce     
      call MPI_REDUCE(xexit,xexitGLOBAL,1,
     +    MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
     
      call MPI_REDUCE(yexit,yexitGLOBAL,1,
     +    MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
     
      call MPI_REDUCE(zexit,zexitGLOBAL,1,
     +    MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
     
      call MPI_REDUCE(deathexit,deathexitGLOBAL,1,
     +    MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
     
      if(id.eq.0)then
      print*,'Average # of scatters per photon:',sngl(nscattGLOBAL
     +                                        /(totphot*numproc))
     
!      print*,'% of fluro photons',real(flucountGLOBAL/
!     +                                nphotons*numproc)*100.

      print*,xexitGLOBAL,yexitGLOBAL,zexitGLOBAL,deathexitGLOBAL

      !write out files
      
      call writer(imageGLOBAL,reflc,trans,Nbins,fileplace,
     +            fluGLOBAL,depositGLOBAL,cbinsnum,fluroglobal
     +            ,cnt,jmeanGLOBAL)
      
      end if
      
      !end MPI processes
      call MPI_Finalize(error)

      end program mcpolar

