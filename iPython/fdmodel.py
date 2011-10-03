from rsf.proj import *

class Dict(dict):
    def __init__(self,**kw):
        dict.__init__(self,**kw)
    
    def __getattr__(self,name):
        if name in self.keys():
            return self[name]
        else:
            return None

    def __setattr__(self,name,val):
        self[name] = val

def defaults(nz,oz,dz,nx,ox,dx,ny=1,oy=0,dy=1,**kw):
    pars = Dict(
        free=False,
        verb=True,
        snap=True,
        expl=False,
        dabc=True,
        cfl=False,
        abcone=False,
        srctype=0,
        fmax=0.0, nb=10, nbell=5,
        jdata=1, jsnap=1, ompnth=4,
        nqz=nz,oqz=oz,dqz=dz,
        nqx=nx,oqx=ox,dqx=dx,
        nqy=ny,oqy=oy,dqy=dy)

    for key,val in kw.items():
        if key in pars.keys():
            pars[key] = val
    return pars

def constant3d(name,val,nz,oz,dz,nx,ox,dx,ny,oy,dy,**kw):
    Flow(name,None,
        '''
        spike n1=%d o1=%f d1=%f
        n2=%d o2=%f d2=%f
        n3=%d o3=%f d3=%f | 
        math output="%f"
        ''' % (nz,oz,dz,nx,ox,dx,ny,oy,dy,val))

def constant2d(name,val,nz,oz,dz,nx,ox,dx,**kw):
    Flow(name,None,
        '''
        spike n1=%d o1=%f d1=%f
        n2=%d o2=%f d2=%f | 
        math output="%f"
        ''' % (nz,oz,dz,nx,ox,dx,val))
        
def point2d(name,z,x,**kw):
    Flow(name,None,
        'spike n1=2 n2=1 nsp=2 k1=1,2 l1=1,2 mag=%f,%f' %(z,x))

def point3d(name,z,x,y,**kw):
    Flow(name,None,
        '''
        spike n1=3 n2=1 nsp=3 k1=1,2,3 l1=1,2,3 mag=%f,%f,%f
        ''' % (z,x,y))

def horizontal2d(name,z,ox,nx,dx,**kw):
    Flow(name+'-x',None,
        '''
        math n1=1 n2=%d o2=%f d2=%f output="x2"
        ''' %(nx,ox,dx))
    Flow(name+'-z',None,
        '''
        math n1=1 n2=%d o2=%f d2=%f output="%f"
        ''' %(nx,ox,dx,z))
    Flow(name,[name+'-z',name+'-x'],
        '''
        cat axis=1 ${SOURCES[1]} 
        ''')

def vertical2d(name,x,oz,nz,dz,**kw):
    Flow(name+'-x',None,
        '''
        math n1=1 n2=%d o2=%f d2=%f output="%f"
        ''' %(nz,oz,dz,x))
    Flow(name+'-z',None,
        '''
        math n1=1 n2=%d o2=%f d2=%f output="x2"
        ''' %(nz,oz,dz))
    Flow(name,[name+'-z',name+'-x'],
        '''
        cat axis=1 ${SOURCES[1]} 
        ''')

def horizontal3d(name,z,ox,nx,dx,oy,ny,dy,**kw):
    Flow(name+'-z',None,
        '''
        math n1=%d n2=%d o1=%f o2=%f d1=%f d2=%f 
        output="%f"
        ''' % (nx,ny,ox,oy,dx,dy,z))
    Flow(name+'-x',name+'-z',
        '''
        math output="x1" | put n1=1 n2=%d
        ''' % (nx*ny))
    Flow(name+'-y',name+'-z',
        '''
        math output="x2" | put n1=1 n2=%d
        ''' % (ny*nx))
    Flow(name,[name+'-'+s for s in ['z','x','y']],
        '''
        put n1=1 n2=%d | 
        cat axis=1 ${SOURCES[1:3]}
        ''' % (nx*ny))

def rec2horizontal3d(name,rec,ox,nx,dx,oy,ny,dy,**kw):
    Flow(name,rec,
        '''
        put n1=%d o1=%f d1=%f n2=%d o2=%f d2=%f
        ''' % (nx,ox,dx,ny,oy,dy))

def waveplot(name,pclip=99,**kw):
    Plot(name,
        '''
        transp plane=12 | 
        graph title="" plotfat=5
        ''')

def ricker(name,freq=20,ot=0,dt=0.001,nt=1001,kt=100,**kw):
    Flow(name,None,
        '''
        spike n1=%(nt)d o1=%(ot)f d1=%(dt)g nsp=1 k1=%(kt)d l1=%(kt)d | 
        ricker1 frequency=%(freq)f | transp plane=12
        ''' % (locals())
        )

def awefd(data,vel,dens,wavelet,sou,rec,wfld=None,
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
           
