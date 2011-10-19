# Utilities for plotting results from elastic data forward
# modeling

from rsf.proj import *


def wavelet(name,pclip=99,**kw):
    Plot(name,
        '''
        transp plane=12 | 
        graph title="" plotfat=5
        ''')

def tensorwavelet(name,pclip=99,**kw):
    Plot(name+'-xx',name,
        '''
        window n2=1 f2=0 |
        graph title="xx" plotfat=5 max2=1.0
        ''')
    Plot(name+'-zz',name,
        '''
        window n2=1 f2=1 |
        graph title="zz" plotfat=5 max2=1.0
        ''')
    Plot(name+'-zx',name,
        '''
        window n2=1 f2=2 |
        graph title="xz" plotfat=5 max2=1.0
        ''')
    Plot(name,[name+'-'+tag for tag in ['xx','zz','zx']],
        'SideBySideAniso')

def data(data,pclip=99,**kw):
    Plot(data+'-z','window n2=1 | grey pclip=%f' % pclip)
    Plot(data+'-x','window f2=1 | grey pclip=%f' % pclip)

def data3d(data,pclip=99,**kw):
    Plot(data+'-z',
        '''
        window n2=1      | 
        grey pclip=%f
        ''' % pclip)
    Plot(data+'-x','window f2=1 n2=1 | grey pclip=%f' % pclip)
    Plot(data+'-y','window f2=2 n2=1 | grey pclip=%f' % pclip)


