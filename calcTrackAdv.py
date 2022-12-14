#!/usr/bin/python3
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt
import numpy as np
from numpy import genfromtxt


data = np.genfromtxt('/home/afoksha/projects/trajectory_calc/trajectory.dat', delimiter=' ')

#ax = plt.axes(projection='3d')


xdata = data[:,0]
ydata = data[:,1]
zdata = data[:,2]
plt.plot(xdata,zdata,'r')
plt.plot(xdata,ydata,'b')
#ax.scatter3D(xdata, ydata, zdata);
plt.title('Projection motion')

plt.show()
