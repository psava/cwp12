from numpy import *

def read(filename,faces=False):

    file = open(filename,'r')
    nvertex, nfacet = fromstring(file.readline(),sep=' ',dtype='int')
    file.close()

    print 'Found %d vertices, %d facets' % (nvertex,nfacet)
    raw_values = loadtxt(filename,skiprows=1)
    vertices = raw_values[:nvertex,1:]
    if faces:
        facets = loadtxt(filename,skiprows=nvertex+1,dtype='int')[:,1:] - 1
        return nvertex, vertices, nfacet, facets
    else:
        return nvertex,vertices

class Arguments(dict):
    def __init__(self,*args):
        dict.__init__(self,args)

    def __setitem__(self,name,val):
        dict.__setitem__(self,name,val)

    def __getitem__(self,name,default=None):
        if name in self.keys():
            return dict.__getitem__(self,name)
        else:
            return default

    def get(self,name,default=None,required=False,output_type=str):
        if name in self.keys():
            return output_type(self[name])
        elif required:
            raise Exception('did not find required argument: %s' % name)
        else:
            return default

    def __str__(self):
        return '\n'.join([key+'='+val for key,val in self.items()])

def initArguments():
    import sys
    args = Arguments()
    for arg in sys.argv[1:]:
#        print 'Getting %s' % arg
        try:
            key,val = arg.split('=')
            args[key] = val
        except Exception, e:
            print e
            print 'ERROR during set'
            pass
    return args

