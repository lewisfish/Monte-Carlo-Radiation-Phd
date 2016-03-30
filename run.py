#!/usr/bin python

import subprocess
import os
import Image

try:
    subprocess.call("./install.sh", shell = True)
    os.chdir(os.getcwd() + '/bin/')
    subprocess.call("./mcgrid", shell = True)

    os.chdir('../data')
#    
#    subprocess.call("gnuplot plot.gp", shell = True)
#    try:
#       Image.open('fluro.png').show()
#    except IOError:
#       print "File does not exist!!"
#       pass
#    os.chdir(os.getcwd() + '/jmean/')
#    subprocess.call('python gui.py', shell=True)
except KeyboardInterrupt:
    pass
