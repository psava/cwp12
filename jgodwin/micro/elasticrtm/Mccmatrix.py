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
    if not '=' in arg and not '.py' in arg:
        print 'found file: ', arg 
        files.append(sf.Input(arg))
    elif 'template' in arg:
        templates.append(sf.Input(arg.split('=')[1]))

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


output = sf.Output()

nfiles = len(files)

output.put('n1',n1)
output.put('n2',n2)
output.put('n3',n3)
output.put('n4',nfiles*nfiles)
output.put('n5',nt)

data = zeros((nfiles,n3,n2,n1),'f')
array = zeros((n3,n2,n1),'f')


for it in range(nt):
    print 'time %04d/%04d' % (it,nt)
    for i in range(nfiles):
        files[i].read(array)
        data[i,:,:,:] = array[:,:,:]

    for j in range(nfiles):
        cc = zeros((nfiles,n3,n2,n1),'f')
        jdata = data[j,:,:,:]
        for k in range(nfiles):
            kdata = data[k,:,:,:]
            maxvals = where(abs(jdata) > abs(kdata), jdata, kdata)
#            cc[k,:,:,:] = (jdata+kdata)**2/(2*(jdata*jdata+kdata*kdata))
            cc[k,:,:,:] = (jdata*kdata)/(maxvals**2)
        output.write(cc)

for file in files:
    file.close()
output.close()
