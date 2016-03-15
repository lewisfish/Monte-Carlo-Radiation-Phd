#!/usr/bin python

import subprocess
import os
import Image

try:
    subprocess.call("./install.sh", shell = True)
    os.chdir(os.getcwd() + '/bin/')
    subprocess.call("mpirun -n 2 ./mcgrid", shell = True)

    os.chdir('../data')
    
        
    subprocess.call("gfortran error.f90 -o error", shell = True)
    subprocess.call("gfortran datareader.f -o data", shell = True)
    subprocess.call("./data")
    subprocess.call("gnuplot plot.gp", shell = True)
    subprocess.call("./error")
    subprocess.call("gnuplot plot.gp", shell = True)
    Image.open('validation.png').show()
    
    os.chdir(os.getcwd() + "/jmean/")
    
    subprocess.call("gfortran error.f90 -o error", shell = True)
    subprocess.call("gnuplot plot.gp", shell = True)
    subprocess.call("./error")
    subprocess.call("gnuplot plot.gp", shell = True)
    Image.open('validation.png').show()
    

    
except KeyboardInterrupt:
    pass
