import numpy
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons
from matplotlib.colors import Normalize
import subprocess
import Tkinter as tk
import tkFileDialog
import sys

root = tk.Tk()
file_path = tkFileDialog.askopenfilename(initialfile='tmp.dat')
root.destroy()


def readslice(inputfilename, ndim):
    shape = (ndim, ndim, ndim, 16)
    fd = open('tmp.dat', 'rb')
    data = numpy.fromfile(
        file=fd, dtype=numpy.float64,sep="").reshape(shape, order='F')
    fd.close()
    data = data[:, :, :, 0]
    return data

print ' '
while True:
    try:
        ndim = raw_input("What is the dimension of the cube?\n")
        break
    except ValueError:
        print 'Please enter an integer!'

subprocess.call("gfortran slicer.f90 -o test", shell=True)
subprocess.call("./test " + file_path + ' ' + ndim, shell=True)

X = readslice('tmp.dat', int(100))

      
class IndexTracker(object):

    def __init__(self, fig, ax, X, view):
        self.ax = ax
        self.view = view
        self.X = X
        self.ax.hold(False)

        rows, cols, self.slices = X.shape
        self.ind = self.slices / 2
        cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.8])
        if(self.view == 'X'):
            self.im = ax.imshow(
                self.X[self.ind, :, :], cmap='cubehelix', interpolation='nearest')
            plt.colorbar(self.im, cax=cbaxes, norm=Normalize(
                vmin=0, vmax=numpy.amax(X[self.ind, :, :])))
        self.update()

    def onclick(self, event):
        if event.button == 1:
            radio.on_clicked(self.viewfunc)

    def viewfunc(self, label):
        self.view = label
        cbaxes = fig.add_axes([0.85, 0.1, 0.03, 0.8])
        if(self.view == 'X'):
            plt.cla()
            self.im = ax.imshow(self.X[self.ind, :, :], cmap='cubehelix')
            plt.colorbar(self.im, cax=cbaxes, norm=Normalize(
                vmin=0, vmax=numpy.amax(X[self.ind, :, :])))
        elif(self.view == 'Y'):
            plt.cla()
            self.im = ax.imshow(self.X[:, self.ind, :], cmap='cubehelix')
            plt.colorbar(self.im, cax=cbaxes, norm=Normalize(
                vmin=0, vmax=numpy.amax(X[:, self.ind, :])))
        elif(self.view == 'Z'):
            self.im = ax.imshow(self.X[:, :, self.ind], cmap='cubehelix')
            plt.colorbar(self.im, cax=cbaxes, norm=Normalize(
                vmin=0, vmax=numpy.amax(X[:, :, self.ind])))
        self.update()

    def onscroll(self, event):
        if event.button == 'up':
            self.ind = numpy.clip(self.ind + 1, 0, self.slices - 1)
        else:
            self.ind = numpy.clip(self.ind - 1, 0, self.slices - 1)
        self.update()

    def update(self):
        if(self.view == 'X'):
            self.im.set_data(self.X[self.ind, :, :])
        elif(self.view == 'Y'):
            self.im.set_data(self.X[:, self.ind, :])
        elif(self.view == 'Z'):
            self.im.set_data(self.X[:, :, self.ind])
        ax.set_title('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()


fig = plt.figure()
ax = fig.add_subplot(111)
fig.canvas.set_window_title('Slices of Fluence')

axcolor = 'lightgoldenrodyellow'
rax = fig.add_axes([0.03, 0.73, 0.15, 0.15],
                   axisbg=axcolor, title='Viewpoint (axis)')

radio = RadioButtons(rax, ('X', 'Y', 'Z'), active=0)

tracker = IndexTracker(fig, ax, X, 'X')

fig.canvas.mpl_connect('scroll_event', tracker.onscroll)
fig.canvas.mpl_connect('button_press_event', tracker.onclick)
plt.show()
