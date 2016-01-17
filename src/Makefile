FCOMP = mpifort
FCFLAGS = -O2
#FCFLAGS = -g -fcheck=all -fbounds-check -fbacktrace 
#-ffpe-trap=zero,overflow,underflow 
#FCFLAGS = -Wall -Wextra -pedantic
#LDFLAGS = 
PROGRAM =  mcgrid 
SRCS =      constants.f95 \
            photon_vars.f95 \
            iarray.f95 \
            opt_prop.f95 \
            subs.f95 \
            ch_opt.f95 \
            density.f95 \
            gridset.f95 \
            noisey.f95 \
            fresnel.f95 \
            taufind1.f95 \
            tauint.f95 \
            stokes.f95 \
            sourceph.f95 \
            peeling.f95 \
            binning.f95 \
            reader.f95 \
            writer.f95 \
            findsmax.f95 \
            mcpolar.f95 \
            ran2.f 
OBJECTS = $(SRCS:.f95=.o)

all:	$(PROGRAM)

$(PROGRAM):	$(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^ 

%.o:  %.f95
	$(FCOMP) $(FCFLAGS) -c $<

.PHONY:	clean

clean:
	rm -f *.o *.mod *.MOD mcgrid
