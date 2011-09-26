import numpy,os,sys

ikeys= ['n%d' % i for i in range(1,10)]
ikeys.append('esize')
fkeys = ['o%d' % i for i in range(1,10)]
fkeys.extend(['d%d' % i for i in range(1,10)])
okeys = ['data_format','in']
okeys.extend(['label%d' % i for i in range(1,10)])
okeys.extend(['unit%d' % i for i in range(1,10)])


class Par(dict):
    def __init__(self):
        dict.__init__(self)
        for arg in sys.argv[1:]:
            try:
                key,val = arg.split('=')
                self[key] = val
            except Exception, e:
                print e
    def bool(self,key,default=None):
        if key in self.keys():
            if self[key] == 'n':
                return False
            elif self[key] == 'y':
                return True
            else:
                return bool(self[key])
        else:
            return default
    def int(self,key,default=None):
        if key in self.keys():
            return int(self[key])
        else:
            return default
    def float(self,key,default=None):
        if key in self.keys():
            return float(self[key])
        else:
            return default
    def string(self,key,default=None):
        if key in self.keys():
            return self[key]
        else:
            return default

class HeaderBinaryError(Exception):
    def __init__(self):
        Exception.__init__(self,'Tried to read a Header/Binary combined file.  Must use only Header/binary separate files with the numpy RSF API')

class Input:
    '''
    An Input is a class to read RSF files.

    Create an Input by passing the name of the rsf file to read, with
    extension attached.

    You may not read header/binary combined files.

    All header parameters are stored in the Input object, and can
    be accessed as attributes of that object.  For example:

    n1 = Input.n1
    n2 = Input.n2
    array = zeros(Input.shape)
    ... and so on
    '''
    def __init__(self,tag='in',*args,**kw):
        if tag == 'in':
            self.header = sys.stdin
        else:
            if tag.endswith('.rsf'):
                name = tag
            else:
                par = Par()
                name = par.string(tag)
            self.header = open(name,'r')

        self.maxdim = 1

        while True:
            try:
                line = self.header.readline()
                if line == '': break
                line = line.strip().strip('\n').strip('\r')
                #print line
            except:
                break
            if '=' in line:
                key,val = line.strip().split('=')
                if key in ikeys:
                    if key.startswith('n') and len(key) == 2:
                        if int(key[1]) > self.maxdim:
                            self.maxdim = int(key[1])
                    setattr(self,key,int(val))
                elif key in fkeys:
                    setattr(self,key,float(val))
                elif key in okeys:
                    setattr(self,key,val.strip('"'))
        if getattr(self,'in') == 'stdin':
            raise HeaderBinaryError
        self.shape = [getattr(self,'n%d' % i,1) for i in range(1,self.maxdim+1)]
        self.shape.reverse()
        if self.shape[0] == 1:
            self.shape.pop(0)
        self.shape = tuple(self.shape)

        #print 'header shape',self.shape

        self.size  = numpy.prod(self.shape)
        
        rdform,rdtype = self.data_format.strip('"').split('_')
        esize = self.esize
        
        if rdtype == 'int':
            self.dtype = 'i%d' % esize
        elif rdtype == 'float':
            self.dtype = 'f%d' % esize
        elif rdtype == 'complex':
            self.dtype = 'c%d' % esize

        self.binary_format = rdform
        if rdform == 'binary':
            binary_mode = 'rb'
        else:
            binary_mode = 'r'
        self.binary = open(getattr(self,'in').strip('"'),binary_mode)

  
    def read(self,shape=None,sep=' '):
        ''' 
        Read the file.  
        
        If shape is specified then only read 
        the number of elements as given in the 
        shape tuple. Return an array with the
        desired shape.
        
        Else, read the full file and return an 
        array with all elements.

        sep is the field delimiter for ASCII files.  
        Default for ASCII files in Madagascar is tab or 
        whitespace delimited.  If you use comma delimited
        then you will need to set this to read properly.
        '''
        
        if shape:
            return self.__read(shape,sep).reshape(shape)
        else:
            return self.__read(self.shape,sep).reshape(self.shape)

    def __read(self,shape,sep=' '):
        '''
        The underlying call to read the actual binary file.

        Do not call yourself, use read() instead.
        '''
        size = numpy.prod(shape)
        if self.binary_format == 'ascii':
            return numpy.fromfile(self.binary,
                dtype=self.dtype,
                count=size,sep=sep)
        else:
            return numpy.fromfile(self.binary,
                dtype=self.dtype,
                count=size)

    def close(self):
        '''
        Close the current Input file now.
        '''
        self.binary.close()
        self.header.close()

    def tell(self):
        '''
        Return the current pointer in the binary file.
        '''
        return self.binary.tell()

    def rewind(self):
        '''
        Reset the binary file to point at the beginning again.
        '''
        self.seek((0,))

    def seek(self,indices):
        '''
        Seek to the given index in the binary file.
        Index should be a tuple, whose indices represent the location
        of the first value that you want to read.  Indices are specified
        in reverse order (i.e. slowest dimension first).  
        
        Remember: Indices
        start from zero and include n-1.
        '''
        if len(indices) > 1:
            element = indices[-1]
            element += numpy.prod(indices[:-1])
        else:
            element = indices[0]
        self.binary.seek(self.esize*element)

    def __str__(self):
        info = ' '.join(['%s=%s' % (key,val) for key,val in self.__dict__.items()])
        return 'Input: %s Closed? %s \n' % (self.header.name,self.header.closed) + info
            
class Output:
    '''
    An Output can be used to write numpy arrays to RSF files.

    You may not write header/binary combined files with this class.

    An Output stores all the header parameters associated for this
    RSF file and then writes them out.  You may set the parameters
    as you would an attribute for any Python object:

    Output.n1 = 12
    Output.d1 = 0.5
    Output.label1='time'
    etc.
    '''
    def __init__(self,tag='out',*args,**kw):
        if tag == 'out':
            self.header = sys.stdout
            stdoutino = os.fstat(sys.stdout.fileno())[1]
            name = None
            for path in os.listdir('.'):
                ino = os.stat(path)[1]
                if ino == stdoutino:
                    name = path
                    break
            if not name: raise Exception('could not find stdout file')
        else:
            if tag.endswith('.rsf'):
                name = tag
            else:
                par = Par()
                name = par.string(tag)
            self.header = open(name,'w')

        path = os.environ.get('DATAPATH')
        setattr(self,'in',os.path.join(path,name+'@'))
        self.binary = open(getattr(self,'in'),'wb')
        self.size = 0 # Number of elements
        self.esize = 0 # Size of elements in bytes
        self.shape = () # Shape of data
        self.nwrites = 0 # Have there been multiple writes?
        self.dtype = None # Type of data numpy.dtype

    def write(self,array):
        '''
        Write the given array to this RSF file.

        Automatically writes the entire array to this file.
        Multiple calls to write may be made, and the
        Output will keep track of how big your last dimension is.

        The first call to write sets the shape of the output file,
        the datatype, and the size of the elements.

        RSF files are regular.  You may only call this method
        with arrays of the same type, and shape.
        '''
        if self.nwrites == 0:
            self.dtype = array.dtype
            self.esize = array.dtype.itemsize
            self.shape = array.shape
            for j in range(1,len(self.shape)+1):
                setattr(self,'n%d' % j, array.shape[-j])

        assert self.shape == array.shape
        assert self.dtype == array.dtype

        array.tofile(self.binary)
        self.size += numpy.prod(array.shape)
        self.nwrites += 1
    
    def getAxesFrom(self,Input):
        '''
        Set the axes to have the same values as the 
        given Input file.
        '''
        for i in range(1,10):
            for pre in ['n','o','d','label','unit']:
                val = getattr(Input,pre+str(i),None)
                if val:
                    setattr(self, pre+str(i), val)
    
    def flush(self):
        '''
        Immediately write any information held in a buffer
        to the header and binary files.

        Typically, you do not need to call this method.  Call close() 
        instead.
        '''
        self.header.flush()
        self.binary.flush()

    def __writeHeader(self):
        '''
        Write the header information given the current state of the
        Output object.  This is called automatically when close() is 
        called.  Do not call this yourself.
        '''
        for key,val in self.__dict__.items():
            if key in ikeys or key in fkeys:
                self.header.write(key+'='+str(val)+'\n')
            elif 'label' in key or 'unit' in key:
                self.header.write(key+'="'+val+'"\n')
        self.header.write('in="'+getattr(self,'in')+'"\n')
        if self.dtype.kind == 'f':
            format = 'float'
        elif self.dtype.kind == 'c':
            format = 'complex'
        else:
            format = 'int'
        self.header.write('data_format="native_' + format +'"\n')
        self.header.write('esize=%s\n' % self.esize)
        if self.nwrites > 1:
            lastdim = len(self.shape)
            self.header.write('n%d=%d\n' % (lastdim+1,self.nwrites))

    def close(self):
        '''
        Close this RSF file, including the header and binary now.

        This is a multi-step process:
            1 - Write the header file.
            2 - Flush both the header and binary buffers.
            3 - Close the Python file identifiers.
        '''
        self.__writeHeader()
        self.flush()
        self.header.close()
        self.binary.close()

    def __str__(self):
        info = ' '.join(['%s=%s' % (key,val) for key,val in self.__dict__.items() if not binary in key])
        return 'Output: %s Closed? %s \n' % (self.header.name,self.header.closed) + info
 
if __name__ == '__main__':

    
    # Test stdin
    input1 = Input()
    data4 = input1.read()
    input1.close()

    # Test stdout
    output5 = Output()
    output5.write(data4)
    output5.getAxesFrom(input1)
    output5.close()

    sys.stdout = sys.__stderr__

    data = numpy.sin(numpy.arange(0,1,.01,'f4'))
    fftd = numpy.fft.fft(data)

    print 'Writing', data.shape
    output = Output('junk.rsf')
    output.label1 = 'time'
    output.unit1  = 'seconds'
    output.write(data)
    output.close()

    output2 = Output('junk2.rsf')
    output2.label1 = 'Frq'
    output2.unit1  = 'Hz'
    output2.write(fftd.astype('c8'))
    output2.close()

    input = Input('junk.rsf')
    print input
    data = input.read()
    print 'mean',data.mean()

    input.seek((50,))
    data2 = input.read((50,))
    print 'mean2',data2.mean(),data2.shape
    input.close()

    output3 = Output('junk3.rsf')
    output3.getAxesFrom(input)
    output3.write(data)
    output3.close()
