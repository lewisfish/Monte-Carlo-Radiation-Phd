#intel
#FCOMP = mpiifort
#FCFLAGS = -O2 -cpp -Dintel
#FCFLAGS = -cpp -Dintel -traceback -check -g
#-gen-interfaces -warn

#gnu
FCOMP = mpifort
FCFLAGS = -O2 -cpp 
#FCFLAGS = -O2 -g -fcheck=all -fbounds-check -fbacktrace 
#-ffpe-trap=zero,overflow,underflow 
#FCFLAGS = -Wall -Wextra -pedantic


#LDFLAGS = 
#goes after $^ -pg 
#-pg goes between fcomp and fcflags
PROGRAM =  mcgrid 

SRCS =      ran2.f \
            constants.f90 \
            photon_vars.f90 \
            iarray.f90 \
            subs.f90 \
            absorbers_mod.f90 \
            opt_prop.f90 \
            fluorophores.f90 \
	    skin_prop.f90 \
            ch_opt.f90 \
            density.f90 \
            gridset.f90 \
            noisey.f90 \
            taufind1.f90 \
            tauint.f90 \
            stokes.f90 \
            sourceph.f90 \
            peeling.f90 \
            binning.f90 \
            reader.f90 \
            writer.f90 \
            findsmax.f90 \
            ppm.f90 \
            mcpolar.f90 
            
OBJECTS = $(SRCS:.f90=.o)

all:	$(PROGRAM)

$(PROGRAM):	$(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^

%.o:  %.f90
	$(FCOMP) $(FCFLAGS) -c $<

.PHONY:	clean

clean:
	rm -f *.o *.mod *.MOD mcgrid
