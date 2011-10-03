from rsf.proj import *
# Plotting utilities for sfacousticfd and sfelasticfd.


def edatagrey2d(data,pclip=99,**kw):
    Plot(data+'-z','window n2=1 | grey pclip=%f' % pclip)
    Plot(data+'-x','window f2=1 | grey pclip=%f' % pclip)

def edatagrey3d(data,pclip=99,**kw):
    Plot(data+'-z',
        '''
        window n2=1      | 
        grey pclip=%f
        ''' % pclip)
    Plot(data+'-x','window f2=1 n2=1 | grey pclip=%f' % pclip)
    Plot(data+'-y','window f2=2 n2=1 | grey pclip=%f' % pclip)

def adatagrey2d(data,pclip=99,**kw):
    pass
