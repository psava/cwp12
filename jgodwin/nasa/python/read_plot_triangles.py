from pylab import *
import sys,triangleio

try:
    name = sys.argv[1]
except:
    raise Exception('Must enter filename to read')

nv,vertices,nf,facets = triangleio.read(name,faces=True)

from enthought.mayavi import mlab
mlab.triangular_mesh(vertices[:,0],vertices[:,1],vertices[:,2],facets)
mlab.show()
