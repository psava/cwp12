from pylab import *
import sys

try:
    filename = sys.argv[1]
except Exception,e:
    raise Exception('Must enter name of file to read.')

file = open(filename,'r')
nvertex, nfacet = fromstring(file.readline(),sep=' ',dtype='int')
file.close()

print 'Found %d vertices, %d facets' % (nvertex,nfacet)
raw_values = loadtxt(filename,skiprows=1)
vertices = raw_values[:nvertex,1:]

facet_map = loadtxt(filename,skiprows=nvertex+1,dtype='int')[:,1:] - 1

from enthought.mayavi import mlab
mlab.triangular_mesh(vertices[:,0],vertices[:,1],vertices[:,2],facet_map)
mlab.show()
