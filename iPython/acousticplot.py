from rsf.proj import *
# Plotting utilities for sfacousticfd and sfelasticfd.


def wavelet(name,pclip=99,**kw):
    Plot(name,
        '''
        transp plane=12 | 
        graph title="" plotfat=5
        ''')

def data(data,pclip=99,**kw):
    pass
