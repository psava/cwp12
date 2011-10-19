class Axis:
    def __init__(self,n,o,d):
        self.n = n
        self.o = o
        self.d = d
        self.min = o
        self.max = o+(n-1)*d

class Model2D:
    '''
    A 2D model is defined with the Z-axis as the fast axis, 
    and the X-axis as the slow axis.  

    Model2D provides a series of commands for creating and
    manipulating data with the same dimensionality.  

    Subclasses define further ideas.
    '''
    def __init__(self,zaxis,xaxis):
        self.z = zaxis
        self.x = xaxis
    
    def homogeneous(self,file,value):
        ''' Create a homogeneous file with the model dimensions  
        and value'''
        Flow(file,None,
            '''
            spike
            n1=%d o1=%f d1=%f
            n2=%d o2=%f d2=%f
            nsp=1 mag=%f
            ''' % (self.z.n,self.z.o,self.z.d,
                self.x.n,self.x.o,self.x.d,value))

    def heterogeneous(self,file,spike_command):
        ''' Create a heterogeneous file using the given command to spike '''
        Flow(file,None,
            '''
            spike
            n1=%d o1=%f d1=%f
            n2=%d o2=%f d2=%f
            ''' % (self.z.n,self.z.o,self.z.d,
                self.x.n,self.x.o,self.x.d,value))

    def view(self,file,**kw):
        ''' View a model file using sfgrey: kw are arguments to pass
        to sfgrey'''
        Result(file,
                '''
                grey 
                parallel2=n labelrot=n wantaxis=y title=""
                pclip=100
                min1=%g max1=%g label1=%s unit1=%s
                min2=%g max2=%g label2=%s unit2=%s
                screenratio=%g screenht=%g wantscalebar=%s
                ''' + 
                ' '.join(map(lambda x: "%s=%s" % (x[0],x[1]),kw)))



