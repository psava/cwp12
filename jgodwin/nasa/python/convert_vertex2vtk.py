from numpy import *
import triangleio,sys

args = triangleio.initArguments()

filename = args.get('in',None,True)
name  = args.get('out',None,True)
nv, vertices = triangleio.read(filename)

file = open(name,'w')
file.write('# vtk DataFile Version 2.0\n')
file.write('Point cloud\n')
file.write('ASCII\n')
file.write('DATASET POLYDATA\n')
file.write('POINTS %d float\n' % nv)
for vertex in vertices:
    file.write('%8f %8f %8f\n' % (vertex[0],vertex[1],vertex[2]))
file.close() 
