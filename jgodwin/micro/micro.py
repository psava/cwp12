from rsf.cluster import *
import random,fdmod


def getpar():

    par = { 
        ###############################
        # Model/Image dimensions
        ###############################
        #'nx':501, 'ox':0, 'dx':0.002,  'lx':'x', 'ux':'km',
        #'ny':151, 'oy':0, 'dy':0.002,  'ly':'y', 'uy':'km',
        #'nz':351, 'oz':0, 'dz':0.002,  'lz':'z', 'uz':'km',
        'nx':251, 'ox':0, 'dx':0.005,  'lx':'x', 'ux':'km',
        'ny':75, 'oy':0, 'dy':0.005,  'ly':'y', 'uy':'km',
        'nz':176, 'oz':0, 'dz':0.005,  'lz':'z', 'uz':'km',
        ###############################
        # Wavelet parameters
        ###############################
        'nt':4001, 'ot':0, 'dt':0.001, 'lt':'t', 'ut':'s',
        'frq':45,      # Peak frequency for Ricker wavelet
        'kt':100,      # Wavelet start position (wavelets are delayd for zero-phase)
        ###############################
        # Modeling code parameters
        ###############################
        'cfl': True,
        'dabc':True,   # Use absorbing boundary condition?
        'nb':80,       # How many cells for absorbing boundary?
        'abcone':True, # Use additional ramp condition for boundaries? (Use default)
        'dsou':False,  # Use displacement source (acoustic-only)
        'expl':False,  # Use exploding reflector (acoustic-only)
        'free':False,  # Use free surface (generate multiples)
        'jdata':1,     # Interval between time-iterations before saving data at recv
        'snap':True,   # Save wavefield snapshots?
        'verb':True,   # Verbose output?
        'jsnap':1,   # Interval between time-iterations before saving wfld snapshot
        'debug':False, # Debugging output (elastic-only)?
        'nbell':5,     # Size of interpolation for injection
        'ssou':False,  # Use stress-source (elastic-only)
        'ompchunk':1,  # OpenMP chunk size (use default)
        'ompnth':2,    # Number of OpenMP threads to use (4 works best)
        ###############################
        # Thomsen parameters for models
        ###############################
        'vp':1.5,
        'ro':2.0,
        ###############################
        # Miscellaneous parameters
        ###############################
        'height':10,
        'nht': 80,
        'nhx': 40,
        'nhz': 40,
        }
    fdmod.param(par)
    par['nframe']=5
    par['iframe']=4
    # ------------------------------------------------------------
    # End user parameters -- NO EDITS BELOW 
    # ------------------------------------------------------------
    par['kz']=2./3.*par['nz']

    return par


def windowreceivers(rr,groups,keys,par):

    for group,gpars in groups.items():
        nwin = gpars['nr']
        owin = gpars['or']
        dwin = gpars['dr']
        Flow('rr-'+group,gpars['group'],'window n2=%d f2=%d j2=%d squeeze=n' % (nwin,owin,dwin))
        Plot('rr-'+group,fdmod.rrplot('plotcol=%d plotfat=10' % gpars['color'],par))

    Flow(rr,['rr-'+group for group in keys],
        'cat axis=2 ${SOURCES[1:%d]}' % len(groups))
    Plot(rr,['rr-'+group for group in keys],'Overlay')


def triangulate(image,tcube,noisy,clean,groups,keys,hypocenters,subgroups,snapshots,par):

    ii = 0
    Fork(nodes=1,time=3,ipn=1)
    for group in keys:
        gpars = groups[group]
        nwin = gpars['nr']
        
        Flow('da-'+group,noisy,
            '''
            window n1=%d f1=%d squeeze=n | 
            ''' % (nwin,ii) + 
            '''put o1=%(oz)g d1=%(dz)g ''' % par)


        Result('da-'+group+'_',clean,
            '''
            window n1=%d f1=%d squeeze=n | 
            ''' % (nwin,ii) + 
            '''
            put o1=0 d1=1 | transp | 
            wiggle poly=y pclip=99 title="" labelsz=6 labelfat=3 titlesz=12 titlefat=3
            label2="\F2 trace\F3 " label1="\F2 time\F3"
            '''   )

        Result('da-'+group,
            '''put o1=0 d1=1 | transp | 
            wiggle poly=n pclip=100 title="" transp=%(transp)d labelsz=6 labelfat=3 titlesz=12 titlefat=3
            yreverse=%(yreverse)d %(custom)s
            label2="\F2 trace\F3 " label1="\F2 time\F3"
            ''' % gpars)

        backproject('da-'+group,'rr-'+group,'vp-2d','ro-2d','_wa-'+group,par)
        Flow('wa-%s'%group,'_wa-%s'%group,
            '''
            transp plane=23 | transp plane=12
            ''' % par)

        Result('wa-%s' % group,'_wa-%s' % group,
            'window f3=%d n3=%d j3=%d | ' % (snapshots[0],snapshots[1],snapshots[2]) + 
            fdmod.cgrey('pclip=100',par))

        for i in range(snapshots[0],snapshots[0]+snapshots[1]*snapshots[2],snapshots[2]):
            Plot('wa-%s-%d' % (group,i),'_wa-%s' % group,'window n3=1 f3=%d | ' % (i) + fdmod.cgrey('pclip=100',par))
            Result('wa-%s-%d' % (group,i),['wa-%s-%d' % (group,i),'rr-2d'],'Overlay')

        subgroupwflds = []
        for sub in subgroups:
            for j in range(0,nwin,sub):
                if sub + j <= nwin:
                    Flow('da-%s-%d-%d' % (group,sub,j),'da-%s' % group,
                        '''
                        window n1=%d f1=%d squeeze=n
                        ''' % (sub,j))
                    Flow('rr-%s-%d-%d' % (group,sub,j),'rr-%s' % group,
                        '''
                        window n2=%d f2=%d squeeze=n
                        ''' % (sub,j))
                else:
                    Flow('da-%s-%d-%d' % (group,sub,j),'da-%s' % group,
                        '''
                        window f1=%d squeeze=n
                        ''' % (j))
                    Flow('rr-%s-%d-%d' % (group,sub,j),'rr-%s' % group,
                        '''
                        window f2=%d squeeze=n
                        ''' % (j))
                backproject('da-%s-%d-%d' % (group,sub,j),
                    'rr-%s-%d-%d' % (group,sub,j),
                    'vp-2d','ro-2d','_wa-%s-%d-%d' % (group,sub,j),par)

                # Go from z-x-t to t-z-x
                Flow('wa-%s-%d-%d'% (group,sub,j),'_wa-%s-%d-%d'%(group,sub,j),
                    '''
                    transp plane=23 | transp plane=12
                    ''' % par)

                Result('wa-%s-%d-%d' % (group,sub,j),'_wa-%s-%d-%d' % (group,sub,j),
                    'window f3=%d n3=%d j3=%d | ' % (snapshots[0],snapshots[1],snapshots[2]) + 
                    fdmod.cgrey('pclip=100',par))

                subgroupwflds.append('_wa-%s-%d-%d' % (group,sub,j))

        j = 0
        for hypocenter in hypocenters:
            xi = hypocenter[0]
            zi = hypocenter[1]
            ti = hypocenter[2]

            Flow('hypo-%d-%s' % (j,group), '_wa-%s' % group,
                '''
                window min1=%(oz)f min2=%(ox)f n1=%(nz)d n2=%(nx)d | 
                ''' % par + 
                '''
                window n1=1 n2=1 f1=%d f2=%d
                ''' % (zi,xi))


            for subgroupwfld in subgroupwflds:
                subgrouphypo = subgroupwfld.replace('_wa','hypo-%d' % j)
                Flow(subgrouphypo,subgroupwfld,
                     '''
                     window n1=1 n2=1 f1=%d f2=%d 
                     ''' % (zi,xi))
            j+= 1


        ii += nwin
        Iterate()

    Join()

    for jhypo in range(len(hypocenters)):
        Flow('hypo-%d' % jhypo, ['hypo-%d-%s' % (jhypo,group) for group in keys],
            '''
            cat axis=2 ${SOURCES[1:%d]}
            ''' % len(keys))
        Result('hypo-%d' % jhypo, 'grey pclip=95')    

        #Save('hypo-%d' % jhypo)
        for sub in subgroups:
            subwflds = []
            for group in keys:
                for j in range(0,groups[group]['nr'],sub):
                    subwflds.append('hypo-%d-%s-%d-%d' % (jhypo,group,sub,j))
            Flow('hypo-%d-%d' % (jhypo,sub), subwflds,
                '''
                cat axis=2 ${SOURCES[1:%d]}
                ''' % len(subwflds))
            Result('hypo-%d-%d' % (jhypo,sub),'grey pclip=95')
            #Save('hypo-%d-%d' % (jhypo,sub))

    for sub in subgroups:
        subwflds = ['wa-%s-%d-%d'% (group,sub,j) for group in keys for j in range(0,groups[group]['nr'],sub) ]

        Flow(tcube+'-sem-%d' % sub,subwflds,
            '''
            semblance m=10 ${SOURCES[1:%d]} | 
            transp plane=12 | transp plane=23
            ''' % len(subwflds))
        Flow(tcube+'-%d' % sub,subwflds,
            '''
            add mode=p ${SOURCES[1:%d]} | 
            transp plane=12 | transp plane=23
            ''' % len(subwflds))
        Result(tcube+'-sem-%d'% sub, 
             'window f3=%d n3=%d j3=%d | ' % 
             (snapshots[0],snapshots[1],snapshots[2]) + 
             fdmod.cgrey('pclip=99.9 gainpanel=a',par))
        Result(tcube+'-%d' % sub, 
            'window f3=%d n3=%d j3=%d | ' % (snapshots[0],snapshots[1],snapshots[2]) + 
            fdmod.cgrey('pclip=99.9 gainpanel=a',par))
        
        Flow(image+'-%d' % sub,tcube+'-%d' % sub,'stack axis=3')
        #Flow(image+'-sem-%d' % sub,tcube+'-sem-%d' % sub,'thr thr=0.4 mode="hard" | stack axis=3')
        Flow(image+'-sem-%d' % sub,tcube+'-sem-%d' % sub,'stack axis=3')
        
        Plot(image+'-sem-box-%d' % sub,image+'-sem-%d' % sub,
            fdmod.cgrey('pclip=100 min2=0.4 max2=0.9 min1=0.2 max1=0.4',par))
        Plot(image+'-box-%d' % sub,image+'-%d' % sub,
            fdmod.cgrey('pclip=99.98 min2=0.4 max2=0.9 min1=0.2 max1=0.4',par))
        Plot(image+'-%d' % sub,fdmod.cgrey('pclip=99.98',par))

        Result(image+'-%d' % sub,[image+'-%d' % sub,'ss-2d','box'],'Overlay')
        Result('image-box-%d' % sub,[image+'-box'+'-%d' % sub,'ss-2d-box'],'Overlay')
        Result('image-sem-box-%d' % sub,[image+'-sem-box-%d' % sub,'ss-2d-box'],'Overlay')
        


    Flow(tcube+'-sem',['wa-%s' % group for group in keys],
            '''
            semblance m=10 ${SOURCES[1:%d]}  | 
            transp plane=12 | transp plane=23
            ''' % len(keys))

    Flow(tcube,['wa-%s'%group for group in keys],
        '''
        add mode=p ${SOURCES[1:%d]} | 
        transp plane=12 | transp plane=23
        ''' % len(keys))
    
    Result(tcube, 
         'window f3=%d n3=%d j3=%d | ' % 
         (snapshots[0],snapshots[1],snapshots[2]) + 
         fdmod.cgrey('pclip=100 gainpanel=a',par))


    Result(tcube+'-sem', 
         'window f3=%d n3=%d j3=%d | ' % 
         (snapshots[0],snapshots[1],snapshots[2]) + 
         fdmod.cgrey('pclip=100 gainpanel=a',par))

    for i in range(snapshots[0],snapshots[0]+snapshots[1]*snapshots[2],snapshots[2]):
        Plot(tcube+'-%d' % i,tcube,'window n3=%d f3=%d | ' % (1,i) + fdmod.cgrey('pclip=99.9 gainpanel=a',par))
        Result(tcube+'-%d' % i , [tcube+'-%d' % i,'rr-2d'],'Overlay')
        Plot(tcube+'-sem-%d' % i,tcube+'-sem','window n3=%d f3=%d | ' % (1,i) + fdmod.cgrey('pclip=99.9 gainpanel=a',par))
        Result(tcube+'-sem-%d' % i , [tcube+'-sem-%d' % i,'rr-2d'],'Overlay')


        
    Flow(image,tcube,'stack axis=3')
    #Flow(image+'-sem',tcube+'-sem','thr thr=0.4 mode="hard" | stack axis=3')
    Flow(image+'-sem',tcube+'-sem','stack axis=3')

    Plot(image+'-box',image,fdmod.cgrey('pclip=99.98 min2=0.4 max2=0.9 min1=0.2 max1=0.4',par))
    Plot(image+'-sem-box',image+'-sem',fdmod.cgrey('pclip=99.98 min2=0.4 max2=0.9 min1=0.2 max1=0.4',par))
    Plot(image,fdmod.cgrey('pclip=99.98',par))
    Plot(image+'-sem',fdmod.cgrey('pclip=100',par))

    Result(image,[image,'ss-2d','box'],'Overlay')
    Result(image+'-sem',[image+'-sem','ss-2d','box'],'Overlay')
    Result('image-box',[image+'-box','ss-2d-box'],'Overlay')
    Result('image-sem-box',[image+'-sem-box','ss-2d-box'],'Overlay')
    

# ------------------------------------------------------------
# Setup functions for calling FD operators
# ------------------------------------------------------------ 
# These operations are usually hidden, but having them here is more 
# transparent.  All possible options are specified by the user.
def backproject(data,receivers,velocity,density,wavefieldname,par):
    Flow(data+'-reversed',data,'sfreverse which=2 opt=i')
    awefd(data+'-junk',wavefieldname,data+'-reversed',
          velocity,density,
          receivers,receivers, par)


def awefd(odat,owfl,idat,velo,dens,sou,rec,par): 
    # call the acoustic wave equation code
    # see sfawe for a more detaile description of options
    
    Flow([odat,owfl],[idat,velo,dens,sou,rec],
         '''
         awe
         ompchunk=%(ompchunk)d ompnth=%(ompnth)d 
         snap=%(snap)d jsnap=%(jsnap)d
         dabc=%(dabc)d nb=%(nb)d
         dsou=%(dsou)d free=%(free)d
         expl=%(expl)d jdata=%(jdata)d
         cfl=%(cfl)d
         fmax=%(frq)f
         verb=%(verb)d
         vel=${SOURCES[1]}
         den=${SOURCES[2]}
         sou=${SOURCES[3]}
         rec=${SOURCES[4]}
         wfl=${TARGETS[1]}
         nqz=%(nz)d
         nqx=%(nx)d
         dqz=%(dz)f
         dqx=%(dx)f
         oqz=%(oz)f
         oqx=%(ox)f
         ''' % par)

# ------------------------------------------------------------
def wavelet(waveletname,frequency,kt,par):
    partemp = par.copy()
    partemp['kt'] = kt
    partemp['frequency'] = frequency
    
    Flow(waveletname,None,
         '''
         spike nsp=1 mag=1 n1=%(nt)d d1=%(dt)g o1=%(ot)g k1=%(kt)d |
         pad end1=%(nt)d |
         ricker1 frequency=%(frequency)g |
         window n1=%(nt)d |
         scale axis=123 |
         put label1=t | thr thr=0.001
         ''' % partemp) 
# ------------------------------------------------------------
def makemicroseisms(ns,wav,sou,par):
    sources = []
    wavelets = []

    r = random.Random()
    r.seed(1234)

    locations = []
    for i in range(ns):
        tag = '-%03d' % i
        xi = r.randrange(100,150)
        zi = r.randrange(50,60)
        ti = r.randrange(par['nt']/4,3*par['nt']/4)

        print 'Microseism %d %d %d %d' % (i,xi,zi,ti)
        locations.append((xi,zi,ti))
        xsou = par['ox']+par['dx']*xi
        zsou = par['oz']+par['dz']*zi
        fdmod.point(sou+tag,xsou,zsou,par)
        wavelet(wav+tag,par['frq'],ti,par)
        sources.append(sou+tag)
        wavelets.append(wav+tag)
        
    Flow(wav+'_',wavelets,'cat axis=2 ${SOURCES[1:%d]}' % ns)
    Flow(sou,sources,'cat axis=2 ${SOURCES[1:%d]}' % ns)

    Plot('ss-2d',fdmod.ssplot('symbol=+ symbolsz=7 plotfat=5',par))
    Plot('ss-2d-box','ss-2d',
        fdmod.ssplot('min1=0.4 max1=0.9 min2=0.2 max2=0.4 plotfat=5 symbol=+ symbolsz=9',par))
    Flow(  'wava','wav_','add scale=10000000 | transp')
    Result('wava','transp |' + fdmod.waveplot('',par))

    # These are bad locations, no microseisms here.
    locations.append((50,25,100))
    locations.append((75,80,100))

    return locations 
# ------------------------------------------------------------
def model(rr,par):
    Flow('zero-2d',None,
         '''
         spike nsp=1 mag=0.0
         n1=%(nz)d o1=%(oz)g d1=%(dz)g 
         n2=%(nx)d o2=%(ox)g d2=%(dx)g |
         put label1=%(lz)s label2=%(lx)s unit1=%(uz)s unit2=%(ux)s
         ''' % par)

    Flow('vz-2d','zero-2d',
        '''
        spike nsp=5
        nsp=5 k1=10,40,70,100,130 l1=39,69,99,129,%(nz)d mag=0.2,0.4,0.6,0.8,1.0 
         n1=%(nz)d o1=%(oz)g d1=%(dz)g 
         n2=%(nx)d o2=%(ox)g d2=%(dx)g |
         put label1=%(lz)s label2=%(lx)s unit1=%(uz)s unit2=%(ux)s | 
        add add=%(vp)f
        ''' % par)

    Flow('fault-2d','zero-2d',
        '''
        spike nsp=1 k1=40 mag=1.0 l1=%(nz)d k2=60 l2=%(nx)d p2=1
         n1=%(nz)d o1=%(oz)g d1=%(dz)g 
         n2=%(nx)d o2=%(ox)g d2=%(dx)g |
         put label1=%(lz)s label2=%(lx)s unit1=%(uz)s unit2=%(ux)s  
        ''' % par)
    Flow('const-2d','zero-2d',
        '''
        spike nsp=1 mag=1.0 k1=40 l1=%(nz)d k2=1 l2=59
         n1=%(nz)d o1=%(oz)g d1=%(dz)g 
         n2=%(nx)d o2=%(ox)g d2=%(dx)g |
         put label1=%(lz)s label2=%(lx)s unit1=%(uz)s unit2=%(ux)s  
        ''' % par)
    Flow('vp-2d','vz-2d','window')

    Flow('ro-2d','zero-2d','math output="%(ro)g"' %par)

    fdmod.makebox('box',0.2,0.4,0.4,0.9,par)
    Plot('box',fdmod.bbplot('',par))

    Plot('vp-2d',fdmod.cgrey('allpos=y pclip=100 bias=1.5 ',par))
    Plot('ro-2d',fdmod.cgrey('bias=2. allpos=y',par))
    Result('vp-2d','vp-2d ss-2d rr-2d box','Overlay')
    Result('ro-2d','ro-2d ss-2d','Overlay')

def synthesize(data,rr,snapshots,par):
    # 2D acoustic modeling
    awefd(data,'wa-2d','wava','vp-2d','ro-2d','ss-2d',rr,par)
    Result(data,'transp |' + fdmod.dgrey('',par))
    for i in range(snapshots[0],snapshots[0]+snapshots[1]*snapshots[2],snapshots[2]):
        Plot('wa-2d-%d' % i,'wa-2d','window n3=%d f3=%d | ' % (1,i) + fdmod.cgrey('pclip=99.9 gainpanel=a',par))
        Result('wa-2d-%d' %i , ['wa-2d-%d' % i,rr],'Overlay')

def addnoise(noisy,data,scale,snapshots,par):
    Flow(noisy,data, 'math output="0" | noise seed=123 | transp | bandpass flo=20 fhi=50 | transp | add scale=%f | add mode=a ${SOURCES[0]} | add scale=1e6' % scale)

    Result(noisy,'transp | grey pclip=99.9')
    backproject(noisy,'rr-2d','vp-2d','ro-2d','wa-%s'% noisy,par)

    Result('wa-%s' % noisy,
        'window f3=%d n3=%d j3=%d | ' % (snapshots[0],snapshots[1],snapshots[2]) + 
        fdmod.cgrey('pclip=100',par))



    
