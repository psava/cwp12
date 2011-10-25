from rsf.proj import *


def defaults(nz,oz,dz,nx,ox,dx,ny=1,oy=0,dy=1,**kw):
    pars = dict(
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

def vofz2d(name,dvdz,v0,nz,oz,dz,nx,ox,dx,**kw):
    Flow(name,None,
        '''
        spike n1=%d o1=%f d1=%f
        n2=%d o2=%f d2=%f | 
        math output="x2*%f+%f"
        ''' % (nz,oz,dz,nx,ox,dx,dvdz,v0))

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
        'spike n1=2 n2=1 nsp=2 k1=1,2 l1=1,2 mag=%f,%f' %(x,z))

def point3d(name,z,x,y,**kw):
    Flow(name,None,
        '''
        spike n1=3 n2=1 nsp=3 k1=1,2,3 l1=1,2,3 mag=%f,%f,%f
        ''' % (x,y,z))

def rec2horizontal2d(output,input,ox,nx,dx,**kw):
    '''
    Put receivers in horizontal grid format again...
    aka change the coordinates from receiver to space.
    '''
    Flow(output,input,
        '''
        put o1=%f d1=%f
        ''' % (ox,dx))

def horizontal2d(name,z,ox,nx,dx,**kw):
    Flow(name+'-x',None,
        '''
        math n1=1 n2=%d o2=%f d2=%f output="x2"
        ''' %(nx,ox,dx))
    Flow(name+'-z',None,
        '''
        math n1=1 n2=%d o2=%f d2=%f output="%f"
        ''' %(nx,ox,dx,z))
    Flow(name,[name+'-x',name+'-z'],
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
    Flow(name,[name+'-x',name+'-z'],
        '''
        cat axis=1 ${SOURCES[1]} 
        ''')

def rec2vertical2d(output,input,oz,dz,**kw):
    Flow(output,input,
        '''
        put o1=%f d1=%f
        ''' % (oz,dz))

def horizontal3d(name,z,ox,nx,dx,oy,ny,dy,**kw):
    Flow(name+'-z_',None,
        '''
        math n1=%d n2=%d o1=%f o2=%f d1=%f d2=%f 
        output="%f"
        ''' % (nx,ny,ox,oy,dx,dy,z))
    Flow(name+'-x',name+'-z_',
        '''
        math output="x1" | put n1=1 n2=%d
        ''' % (nx*ny))
    Flow(name+'-y',name+'-z_',
        '''
        math output="x2" | put n1=1 n2=%d
        ''' % (ny*nx))
    Flow(name+'-z',name+'-z_',
        '''
        put n1=1 n2=%d
        ''' % (nx*ny))
    Flow(name,[name+'-'+s for s in ['x','y','z']],
        '''
        cat axis=1 ${SOURCES[1:3]}
        ''' )

def rec2horizontal3d(name,rec,ox,nx,dx,oy,ny,dy,**kw):
    Flow(name,rec,
        '''
        put n1=%d o1=%f d1=%f n2=%d o2=%f d2=%f
        ''' % (nx,ox,dx,ny,oy,dy))


