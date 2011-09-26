from numpy import *
import sys, triangleio

args = triangleio.initArguments()

name = args.get('in',None,required=True)
outname = args.get('out',None,required=True)
tol  = args.get('tol',0.1,False,float)
maxpoints = args.get('maxpts',100000,False,int)

f = open(name,'r')
for i in range(0,5):
    line = f.readline()
n1,n2,n3 = [int(x) for x in line.split()[1:]]

for i in range(5,10):
    line = f.readline()
line = f.read()
f.close()

print 'Dimensions: %d %d %d ' % (n1,n2,n3)
data = fromstring(line,sep=' ')
print 'Shape check: ', data.shape[0] == n1*n2*n3

data[where(data > tol)] = -1
data[where(data < -tol)] = -1
data[where(data > -tol)] = 1
data[where(data < 0 )] = 0

data = data.reshape((n3,n2,n1))

seed = (n3/2,n2/2,n1/2)

iteration = 0; maxiter = 10000
converged = False
npoints = 1
boundary_points = zeros((3,maxpoints),dtype='i')
boundary_points[:,0] = n3/2,n2/2,n1/2

while not converged and iteration < maxiter:
    print "\rIteration: %05d %05d" % (iteration,npoints)
    temp_points     = zeros((3,maxpoints),dtype='i')
    for ip in range(npoints):
        i3,i2,i1 = boundary_points[:,ip]
        data[i3,i2,i1] = 1
    new_points = 0
    for ip in range(npoints):
        i3,i2,i1 = boundary_points[:,ip]
        for ii3 in range(i3-1,i3+2,1):
            if ii3 < n3 and ii3 >=0:
                for ii2 in range(i2-1,i2+2,1):
                    if ii2 < n2 and ii2 >=0:
                        for ii1 in range(i1-1,i1+2,1):
                            if ii1 < n1 and ii1 >=0 and data[ii3,ii2,ii1] == 0:
                                data[ii3,ii2,ii1] = 1
                                temp_points[:,new_points] = ii3,ii2,ii1
                                new_points += 1

    npoints = new_points
    if npoints == 0:
        break
    boundary_points = temp_points[:,:]
    iteration += 1

data = data.astype('int').reshape((n3*n2,n1))
savetxt(outname,data,fmt='%03d',delimiter=' ')
