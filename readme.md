#        3D SCATTERED LIGHT CODE

This directory contains the source codes for the 3D Cartesian grid.
At the end of the simulation the code outputs to the screen the average number of scatterings 
per Monte Carlo photon packet, and fluences into the data folder.

#### The FORTRAN source files are:

            density.f95 
            gridset.f95 
            iarray.f95 
            mcpolar.f95
            ran2.f 
            findsmax.f95 
            stokes.f95 
            sourceph.f95
            noisey.f95 
            fresnel.f95
            tauint2.f95
            peeling.f95 
            taufind1.f95
            binning.f95 
            reader.f95
            writer.f95

#### The FORTRAN modules are:

            constants.mod 
            photon.mod
            subs.mod

#### Input parameters are in:

	input.params
	opt.params
	noisedots.dat

##### The file that compiles the code and creates the executable file 'mcgrid' is:

	install.sh
	
	this can be run by ./install.sh
	may have to change permissions first in order to execute the script.
	This can be done by using sudo chmod +755 instsll.sh on linux


