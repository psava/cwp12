#!/usr/bin/env python
import rsf.api as sf
from numpy import *
import sys

sys.stdout = sys.__stderr__

par = sf.Par()

#input = sf.Input()

files = []

templates = []
for arg in sys.argv:
    if 'template' in arg:
        print 'found template:', arg
        templates.append(sf.Input(arg))
    elif not '=' in arg and not '.py' in arg:
        print 'found file: ', arg 
        files.append(sf.Input(arg))

print 'found %d files and %d templates' % (len(files),len(templates))
n1 = files[0].int('n1')
n2 = files[0].int('n2')
n3 = files[0].int('n3')
n4 = files[0].int('n4')

if n4 == 1:
    nt = n3
    dim = 2
    n3 = 1
else:
    nt = n4
    dim = 3
    n3 = n3

nfiles = len(files)
matrix = par.string('matrix',None)

image  = sf.Output()
image.put('n1',n1)
image.put('n2',n2)
image.put('n3',n3)
image.put('n4',len(templates))
image.put('n5',nt)


if matrix:
    output = sf.Output('matrix')
    output.put('n1',n1)
    output.put('n2',n2)
    output.put('n3',n3)
    output.put('n4',nfiles*nfiles)
    output.put('n5',nt)

data = zeros((nfiles,n3,n2,n1),'f')
array = zeros((n3,n2,n1),'f')

ntemp = len(templates)
templatearray = zeros((ntemp,nfiles,nfiles),'f')
for itemp in range(ntemp):
    template = templates[itemp]
    template.read(templatearray[itemp,:,:])

print 'read templates', templatearray.min(), templatearray.max()

for it in range(nt):
    print 'time %04d/%04d' % (it,nt)
    for i in range(nfiles):
        files[i].read(array)
        data[i,:,:,:] = array[:,:,:]

    cc = zeros((nfiles*nfiles,n3,n2,n1),'f')
    for j in range(nfiles):
        jdata = data[j,:,:,:]
        for k in range(nfiles):
            ii = j*nfiles+k
            kdata = data[k,:,:,:]
            maxvals = where(abs(jdata) > abs(kdata), jdata, kdata)
#            cc[k,:,:,:] = (jdata+kdata)**2/(2*(jdata*jdata+kdata*kdata))
            cc[ii,:,:,:] = (jdata*kdata)/(maxvals**2)
            #cc[ii,:,:,:] = sign(cc[ii,:,:,:])
            cc[ii,:,:,:] = cc[ii,:,:,:]

    if matrix:
        output.write(cc)

    print 'computing template scores'
    for itemp in range(ntemp):
        score = zeros((n3,n2,n1),'f')
        for j in range(nfiles):
            for k in range(nfiles):
                tempval = templatearray[itemp,j,k]
                score[:,:,:] += tempval*cc[j*nfiles+k,:,:,:]
        image.write(score)
    
for file in files:
    file.close()
for template in templates:
    template.close()
if matrix:
    output.close()
image.close()
