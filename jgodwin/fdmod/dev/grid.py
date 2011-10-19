from pylab import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

np = 4
X,Y,Z = mgrid[0:(np-1):1j*np,0:(np-1):1j*np,0:(np-1):1j*np]

ax.scatter(X.flatten(),Y.flatten(),Z.flatten())


sp = np - 1

X1,Y1,Z1 = mgrid[0.5:sp-0.5:1j*sp,0.5:sp-0.5:1j*sp,0.5:sp-0.5:1j*sp]

ax.scatter(X1.flatten(),Y1.flatten(),Z1.flatten(),c='r',marker='^')



show()

