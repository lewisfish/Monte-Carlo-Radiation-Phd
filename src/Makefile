FCOMP = mpifort
FCFLAGS = -Ofast
#FCFLAGS = -g -fcheck=all -fbounds-check -fbacktrace 
#-ffpe-trap=zero,overflow,underflow 
#FCFLAGS = -Wall -Wextra -pedantic
#LDFLAGS = 
PROGRAM =  mcgrid 
SRCS =      constants.f95 \
            photon.f95 \
            subs.f95 \
            density.f95 \
            gridset.f95 \
            iarray.f95 \
            mcpolar.f95 \
            ran2.f \
            findsmax.f95 \
            stokes.f95 \
            sourceph.f95 \
            noisey.f95 \
            fresnel.f95 \
            tauint2.f95 \
            peeling.f95 \
            taufind1.f95 \
            binning.f95 \
            reader.f95 \
            writer.f95
OBJECTS = $(SRCS:.f95=.o)

all:	$(PROGRAM)

$(PROGRAM):	$(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^ 

%.o:  %.f95
	$(FCOMP) $(FCFLAGS) -c $<

.PHONY:	clean

clean:
	rm -f *.o *.mod *.MOD mcgrid
