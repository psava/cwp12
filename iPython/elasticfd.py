from rsf.proj import *


def ricker(name,freq=20,ot=0,dt=0.001,nt=1001,kt=100,**kw):
    # An explosive source in 2D
    Flow(name,None,
        '''
        spike n1=%(nt)d o1=%(ot)f d1=%(dt)g nsp=1 k1=%(kt)d l1=%(kt)d | 
        ricker1 frequency=%(freq)f | transp plane=12 | spray axis=2 n=2
        ''' % (locals())
        )

def ricker3d(name,freq=20,ot=0,dt=0.001,nt=1001,kt=100,**kw):
    # An explosive source in 3D
    Flow(name,None,
        '''
        spike n1=%(nt)d o1=%(ot)f d1=%(dt)g nsp=1 k1=%(kt)d l1=%(kt)d | 
        ricker1 frequency=%(freq)f | transp plane=12 | spray axis=2 n=3
        ''' % (locals())
    )

def rickerdc2(name,freq=20,ot=0,dt=0.001,nt=1001,kt=100,**kw):
    ricker(name+'_',freq,ot,dt,nt,kt)
    Flow(name+'_mask',None,
        '''
        spike n1=%(nt)d d1=%(dt)f n2=2 k2=1,2 l2=1,2 nsp=2 mag=1.0,-1.0 | 
        transp plane=12 | transp plane=23 | transp plane=12
        ''' % locals())
    Flow(name,[name+'_',name+'_mask'],
        '''
        add mode=p ${SOURCES[1]}
        ''')

def rickershear(name,mu=1.0,D=0.1,S=0.01,freq=20,ot=0,dt=0.001,nt=1001,kt=100,**kw):
    Flow(name,None,
        '''
        spike n1=%(nt)d o1=%(ot)f d1=%(dt)g nsp=1 k1=%(kt)d l1=%(kt)d | 
        ricker1 frequency=%(freq)f | transp plane=12 | spray axis=2 n=3 | 
        cut n2=2 
        ''' % locals())


def rickerdc(name,theta=0.0,mu=1.0,D=0.1,S=.01,
    freq=20,ot=0,dt=0.001,nt=1001,kt=100,**kw):
    # A Double couple ricker wavelet, use with tensor source
    # type in the elastic modeling code (srctype=3) 
    # Theta is the rotation of the fault, WRT to the horizontal
    Flow(name+'_',None,
        '''
        spike n1=%(nt)d o1=%(ot)f d1=%(dt)g nsp=1 k1=%(kt)d l1=%(kt)d | 
        ricker1 frequency=%(freq)f  |  transp plane=12 | transp plane=23 | 
        scale axis=123
        ''' % locals())
    import math
    rads = math.radians(theta)

    xxscale = 2*math.sin(rads)*math.cos(rads)
    zxscale = math.sin(rads)**2+math.cos(rads)**2

    Flow(name+'_xx',name+'_',
        '''
        math output="input*%f"
        ''' % (xxscale))

    Flow(name+'_zz',name+'_',
        '''
        math output="input*%f"
        ''' % (xxscale))

    Flow(name+'_zx',name+'_',
        '''
        math output="input*%f"
        ''' % (zxscale))

    Flow(name,[name+'_'+tag for tag in ['xx','zz','zx']],
        '''
        cat axis=2 ${SOURCES[1:3]}
        ''')

def windowDensity(name,dens,**kw):
    Flow(name,dens,
        '''
        window f1=1 f2=1 f3=1
        ''')

def run(data,ccc,dens,wavelet,sou,rec,wfld=None,
    verb=False,snap=False,free=False, expl=False, dabc=True, 
    cfl=False, abcone=False, srctype=0, ani=None,
    fmax=0.0,jdata=1,jsnap=1,
    nb=10, ompnth=8, nbell=5,
    nqz=1,oqz=0,dqz=1,
    nqx=1,oqx=0,dqx=1,
    nqy=1,oqy=0,dqy=1,**kw):

    if ani == None:
        raise Exception('Must specify the type of anisotropy')

    Flow([data,wfld],[wavelet,ccc,dens,sou,rec],
        '''
        /opt/rsf/bin/sfewe_fixed_dev
        wfl=${TARGETS[1]}
        ccc=${SOURCES[1]}
        den=${SOURCES[2]}
        sou=${SOURCES[3]}
        rec=${SOURCES[4]}
        ani=%(ani)d
        ompnth=%(ompnth)d nbell=%(nbell)d
        verb=%(verb)d free=%(free)d
        expl=%(expl)d dabc=%(dabc)d nb=%(nb)d
        cfl=%(cfl)d fmax=%(fmax)f 
        abcone=%(abcone)d srctype=%(srctype)d 
        snap=%(snap)d jdata=%(jdata)d jsnap=%(jsnap)d
        nqx=%(nqx)d nqz=%(nqz)d nqy=%(nqy)d
        oqz=%(oqz)f oqx=%(oqx)f oqy=%(oqy)f
        dqy=%(dqy)f dqz=%(dqz)f dqx=%(dqx)f
        '''  % (locals()))


