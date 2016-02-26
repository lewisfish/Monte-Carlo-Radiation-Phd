#        3D SCATTERED LIGHT MONTE CARLO CODE

This repository contains the source codes for the 3D grid Cartesian monte carlo radiation transfer codes developed as part of my phd.
The codes simulate the transfer of light through tissue of various types and optical properties.
Currently there are two codes, the master and the spain branches. The spain branch has a cuboid of Nd:yag under tissue that fluoresces at three diffrent wavelengths.
The master branch is more general and can be adapted for purpose.
At the end of the simulation the code outputs to the screen the average number of scatterings 
per Monte Carlo photon packet, and fluences into the data folder.

#### The FORTRAN source files are:

            constants.f95        contains the various constants used in the simulation
            photon_vars.f95      contains the various photons variables
            iarray.f95           contains the arrays variables names
            subs.f95             inits arrays and directory paths 
            density.f95          sets the grids density
            noisey.f95           calculates the 'real' surface normal due to bumpmapping
            fresnel.f95          calculates fresnel coeff.
            taufind1.f95         calculates tau in a direction
            tauint.f95           generates a tau value and integrates through the grid
            gridset.f95          sets up the grid
            stokes.f95           scatters the photon
            sourceph.f95         photon source subroutine
            peeling.f95          peeling off subroutine
            binning.f95          binning subroutine
            reader.f95           reads in optical parameters
            writer.f95           writes out results
            findsmax.f95         finds smax need for tauint calculation
            mcpolar.f95          main program
            ran2.f               random number generator

#### Input parameters are in:

	input.params
	opt.params
	noisedots.dat

##### The file that compiles the code and creates the executable file 'mcgrid' is:

   Only been tested on linux so far. works with intel and gfortran compilers.
   Also been tested on computing clusters. http://www-solar.mcs.st-and.ac.uk/~herbert/cluster/
	
	install.sh
	
	This can be run by ./install.sh
	May have to change permissions first in order to execute the script.
	This can be done by using sudo chmod +755 instsll.sh on linux

#### To do/implement
   finish/test scatt_method and fix binning to be more general
   sort fresnel so can differ between transmission on diffrent faces not just top and bottom/ tidy code up
   investigate +1 in celli,j,k calc in tauint2
   fix up fluro stuff i.e. give options for what type of fluro being done
   change opt arrays to allocatble (as part of above point)
   test bump map, set cross as noise shape then look at diffrent slices/ do calculations manualy
   add proper treatment of peeling off for noisy surface
   add colours/render as option
   proper formatting
   parallize-done but not 100% happy with. change makefile so that mpi is an option.
