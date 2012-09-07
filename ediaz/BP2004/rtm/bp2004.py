try:
    from rsf.cluster import *
except:
    from rsf.proj import *
import fdmod

# ------------------------------------------------------------
# model parameters
def param():
    par = {
        'nx':5395,  'ox':0.000, 'dx':0.0125,  'lx':'x',       'ux':'km',
        'ny':1,     'oy':0.000, 'dy':0.01250,  'ly':'y',       'uy':'km',
        'nz':1911,  'oz':0,     'dz':0.00625,  'lz':'z',       'uz':'km',
        'nt':2001,  'ot':0,     'dt':0.006,    'lt':'t',       'ut':'s',
        'nh':1201,  'oh':0.0,   'dh':0.0125,   'lh':'offset', 'uh':'km',
        'ns':1348 , 'os':1.0,   'ds':0.050 ,   'ls':'offset', 'us':'km'
        }
    
    par['nb']=250


    return par

# ------------------------------------------------------------
def modpar(par):

    par['kt']=100
    par['nt']=12001
    par['dt']=0.001
    par['nb']=150
    par['jsnap']=1000
    par['jdata']=1
    par['wweight']=50
    par['wclip']=0.5
 
# ------------------------------------------------------------   
def migpar(par):

    par['verb']='y'
    par['eps']=0.1
    par['nrmax']=5
    par['dtmax']=0.00005
    par['tmx']=16
    
    par['fw']=36
    par['jw']=1
    par['dw']=1/(par['nt']*par['dt'])
    par['kw']=par['nt']/2+1
    par['ow']=par['fw']*par['dw']
    par['nw']=240
    par['eic']='itype=o'

# ------------------------------------------------------------
def get_vellw(velo,par):

    migvelfile = 'data/bp2004/vel_z6.25m_x12.5m_lw.segy'
    #Fetch(velo,'sigsbee')

    Flow([velo+'-raw',velo+'-t','./'+velo+'-h','./'+velo+'-b'],
         migvelfile,
         '''
         segyread
         tape=$SOURCE
         tfile=${TARGETS[1]}
         hfile=${TARGETS[2]}
         bfile=${TARGETS[3]}
         ''',stdin=0)
    
    Flow(velo,
         velo+'-raw',
         '''
         scale rscale=0.001 |
         put
         o1=%(oz)g d1=%(dz)g label1=%(lz)s unit1=%(uz)s
         o2=%(ox)g d2=%(dx)g label2=%(lx)s unit2=%(ux)s
         ''' %par )
# ------------------------------------------------------------
def get_velexact(velo,par):

    migvelfile = 'data/bp2004/vel_z6.25m_x12.5m_exact.segy'
    #Fetch(velo,'sigsbee')

    Flow([velo+'-raw',velo+'-t','./'+velo+'-h','./'+velo+'-b'],
         migvelfile,
         '''
         segyread
         tape=$SOURCE
         tfile=${TARGETS[1]}
         hfile=${TARGETS[2]}
         bfile=${TARGETS[3]}
         ''',stdin=0)
    
    Flow(velo,
         velo+'-raw',
         '''
         scale rscale=0.001 |
         put
         o1=%(oz)g d1=%(dz)g label1=%(lz)s unit1=%(uz)s
         o2=%(ox)g d2=%(dx)g label2=%(lx)s unit2=%(ux)s
         ''' %par )
# ------------------------------------------------------------
def get_velexact_ns(velo,par):

    migvelfile = 'data/bp2004/vel_z6.25m_x12.5m_nosalt.segy'
    #Fetch(velo,'sigsbee')

    Flow([velo],
         migvelfile,
         '''
         segyread
         tape=$SOURCE | 
         scale rscale=0.001 |
         put
         o1=%(oz)g d1=%(dz)g label1=%(lz)s unit1=%(uz)s
         o2=%(ox)g d2=%(dx)g label2=%(lx)s unit2=%(ux)s
         ''' %par,stdin=0 )

# ------------------------------------------------------------
def get_denexact(den,par):

    denfile = 'data/bp2004/density_z6.25m_x12.5m.segy'
    #Fetch(velo,'sigsbee')

    Flow([den],
         denfile,
         '''
         segyread
         tape=$SOURCE | 
         put
         o1=%(oz)g d1=%(dz)g label1=%(lz)s unit1=%(uz)s
         o2=%(ox)g d2=%(dx)g label2=%(lx)s unit2=%(ux)s 
         ''' %par,stdin=0 )



# ------------------------------------------------------------
def make_shots_0001_0200(shots,par):
    data='data/bp2004/shots0001_0200.segy'
    
    Flow(shots,data,
        '''
        segyread tape=${SOURCES[0]} |
        put o3=%(os)g d3=%(ds)g n3=200
                                n2=%(nh)g 
            o1=%(ot)g d1=%(dt)g n1=%(nt)d|
        reverse opt=i which=2 |
        put o2=%(oh)g d2=%(dh)g n2=%(nh)g label2=%(lh)s 
        '''%par,stdin=0)


# ------------------------------------------------------------
def make_shots_1201_1348(shots,par):
    data='data/bp2004/shots1201_1348.segy'
    
    Flow(shots,data,
        '''
        segyread tape=${SOURCES[0]} |
        put o3=%(os)g d3=%(ds)g n3=148
                                n2=%(nh)g 
            o1=%(ot)g d1=%(dt)g n1=%(nt)d|
        reverse opt=i which=2 |
        put o2=%(oh)g d2=%(dh)g n2=%(nh)g label2=%(lh)s |
        '''%par+
        '''
        put o3=%g
        '''%((1200)*par['ds']+par['os']),stdin=0)



# ------------------------------------------------------------
def make_rr(rr,shot,par):
    
    xs=(shot-1)*par['ds']+par['os']

    Flow(rr,None,
        '''
        spike n1=%(nh)d d1=%(dh)g o1=%(oh)g | 
        rtoc | '''%par+'''
        math output="real(%g-x1) +imag(0)"|
        dd type=float |'''%xs+'''
        put n1=2 n2=%(nh)d 
        '''%par)

# ------------------------------------------------------------
def make_ss(ss,shot,par):
    
    xs=(shot-1)*par['ds']+par['os']

    Flow(ss,None,
        '''
        spike n1=1 | 
        rtoc | '''%par+'''
        math output="real(%g) +imag(0)"|
        dd type=float |'''%xs+'''
        put n1=2 n2=1 
        '''%par)

