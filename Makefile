#*************** Makefile created by Mike Wolff ****************************

#******************************** G77/Linux Fortran ************************
FC     =       gfortran -fbounds-check
#EXTRA_OPT =     -mpentium -malign-double -fforce-mem -fforce-addr \
#                -ffast-math -funroll-all-loops
# May want to experiment by adding the extra optimization flags to get
## better runtime. But then again, maybe not.
#FFLAGS  =       -O2 $(EXTRA_OPT) -ffloat-store
FFLAGS  =       -Ofast
LDFLAGS = 
time_it         = get_cpu_sun

#******************************** PGI Fortran ************************
#FC      =       pgf77
#FFLAGS  =      -fast
#LDFLAGS =	-fast 
#time_it         = get_cpu_sun

#******************************** Sun Fortran ************************
#FC     =       f77
#FFLAGS  =      -fast -O
#LDFLAGS =	-fast -O
#time_it         = get_cpu_sun

#******************************** Lahey-Fujitsu lf95 ************************
#
#FC      =       lf95
#FFLAGS  =       --tpp --nsav -O --nwarn -c
#LDFLAGS =
#time_it         = get_cpu_sun

#****************************************************************************


OBJSB	=	density.o \
	     	gridset.o \
            iarray.o \
            mcpolar.o \
            ran2.o \
            stokes.o \
            sourceph.o \
            noisey.o \
            fresnel.o \
            tauint2.o \
            peeling.o \
            gridedge.o \
            taufind1.o \
            writer.f

mcgrid:	$(OBJSB)
		$(FC) $(OBJSB) $(LDFLAGS) -o mcgrid

clean:;		/bin/rm -f *.o

