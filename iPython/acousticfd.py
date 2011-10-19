from rsf.proj import *

def wfldgrey(wfld,pclip=100,gainpanel='a',result=False,**kw):
    if result:
        Result(wfld,
            '''
            ''')
    else:
        Plot(wfld,
            '''
            ''')

def ricker(name,freq=20,ot=0,dt=0.001,nt=1001,kt=100,**kw):
    Flow(name,None,
        '''
        spike n1=%(nt)d o1=%(ot)f d1=%(dt)g nsp=1 k1=%(kt)d l1=%(kt)d | 
        ricker1 frequency=%(freq)f | transp plane=12
        ''' % (locals())
        )

def run(data,vel,dens,wavelet,sou,rec,wfld=None,
    verb=False,snap=False,free=False, expl=False, dabc=True, 
    cfl=False, abcone=False, srctype=0, fmax=0.0,jdata=1,jsnap=1,
    nb=10, ompnth=8, nbell=5,
    nqz=1,oqz=0,dqz=1,
    nqx=1,oqx=0,dqx=1,
    nqy=1,oqy=0,dqy=1,time=False,binary=None,**kw):

    targets = [data]
    if wfld:
        targets.append(wfld)
    sources = [wavelet,vel,dens,sou,rec]

    print locals()
    if not time:
        Flow(targets,sources,
            '''
            awe
            wfl=${TARGETS[1]}
            vel=${SOURCES[1]}
            den=${SOURCES[2]}
            sou=${SOURCES[3]}
            rec=${SOURCES[4]}
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
    else:
        time_bin = WhereIs('time')
        if not binary: binary = 'sfawe'
        binary = WhereIs(binary)
        Flow(targets,sources,
            '''
            %(time_bin)s 
            %(binary)s
            wfl=${TARGETS[1]}
            vel=${SOURCES[1]}
            den=${SOURCES[2]}
            sou=${SOURCES[3]}
            rec=${SOURCES[4]}
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


