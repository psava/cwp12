from pylab import *

'''
From Calculation of weights in finite difference formulas.
Bengt Fornberg
'''
scale = 1.0
m  = 1
z = 0
nd = 6
n = nd-1

dx = 1.0
x  = arange(-n/2+0.5,n/2+1.5,1) # Use half-distance for ewe

c  = zeros((nd,m+1))
c1 = 1.0
c4 = x[0]-z

c[0,0] = 1.0
for i in range(1,n+1):
    mn = min(i,m)
    c2 = 1.0
    c5 = c4
    c4 = x[i] - z
    for j in range(0,i):
        c3 = x[i] - x[j]
        c2 = c2*c3
        if j == i-1: 
            for k in range(mn,0,-1):
                c[i,k] = c1*(k*c[i-1,k-1]-c5*c[i-1,k])/c2
            c[i,0] = -c1*c5*c[i-1,0]/c2
        for k in range(mn, 0, -1):
            c[j,k] = (c4*c[j,k] - k*c[j,k-1])/c3
        c[j,0] = c4*c[j,0]/c3
    c1 = c2

for i in range(c.shape[1]):
    print '--- Derivative %02d ---' % i
    for j in range(c.shape[0]):
        print 'C%d - %08f' % (j-n/2,scale*c[j,i])
