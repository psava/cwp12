#!/usr/env/bin python
import os,platform,sys,glob,urllib,subprocess,shutil


i586 = 'http://download.java.net/media/jogl/builds/archive/jsr-231-1.1.1a/jogl-1.1.1a-linux-i586.zip'
amd64 = 'http://download.java.net/media/jogl/builds/archive/jsr-231-1.1.1a/jogl-1.1.1a-linux-amd64.zip'
mac = 'http://download.java.net/media/jogl/builds/archive/jsr-231-1.1.1a/jogl-1.1.1a-macosx-universal.zip'

tarch = platform.architecture()[0]
plat = sys.platform
if plat == 'darwin':
    guess = '3, mac'
elif 'linux' in plat:
    if '64' in tarch:
        guess = '2, amd64'
    else:
        guess = '1, i586'
else:
    guess = "* DO NOT KNOW *"
print \
'''
What is your machine architecture?

Choices are:  
    1 i586 (32-bit Linux OS)
    2 amd64 (64-bit Linux OS)
    3 mac  (All Macs)

Hint: Python thinks you should use: %s
''' % guess
arch = raw_input()
choices = {'1':i586,'2':amd64,'3':mac}


cwd = os.getcwd()

try:
    url = choices[arch]
except KeyError, e:
    print 'You did not enter a valid choice!'
    sys.exit(1)

print '-------------------'
print '...Getting zip file from internet'
urllib.urlretrieve(url,'temp_jogl.zip')
print '...Extracting zip file '

p = subprocess.Popen('unzip temp_jogl.zip',shell=True)
p.wait()


jogldir = glob.glob('jogl*')[0]
shutil.move(jogldir,'jogl')

print '...Moving libraries and JARs'
jars = glob.glob('jogl/lib/*.jar')
for jar in jars:
    shutil.move(jar,'jar')


if not os.path.exists('lib'):
    os.mkdir('lib')

if arch == '3':
    libs = glob.glob('jogl/lib/*.jnilib')
else:
    libs = glob.glob('jogl/lib/*.so')

for lib in libs:
    shutil.move(lib,'lib')

print '...Deleting temporary files and folders'
os.remove('temp_jogl.zip')
shutil.rmtree('jogl')
print '...Completed installation'


javahome = cwd

rsf = os.environ.get('RSFROOT')
rsfjar = os.path.join(rsf,'lib','rsf.jar')

if arch == '3':
        lib = '$JAVA_SDK/lib'
        bindir = '$JAVA_SDK/bin'
        jars = map(lambda x: '$JAVA_SDK/%s' % x, glob.glob('jar/*.jar'))
        JAVAHOME='JAVA_SDK'
        LDLIBRARY='DYLD_LIBRARY_PATH'
else:
        lib = '$JAVA_SDK/lib'
        bindir = '$JAVA_SDK/bin'
        jars = map(lambda x: '$JAVA_SDK/%s' % x, glob.glob('jar/*.jar'))
        JAVAHOME='JAVA_SDK'
        LDLIBRARY='LD_LIBRARY_PATH'

jarpath = ':'.join(jars)
classpath = ':'.join([jarpath,bindir,rsfjar,'.'])
ldpath = ':'.join([os.path.join(rsf,'lib'),lib])

print \
''' 
Please set the following environmental variables as follows (COPY/PASTE):

export %s=%s
export CLASSPATH=%s
export %s=%s

PLEASE RECONFIGURE, AND REINSTALL MADAGASCAR USING: API=java option
''' % (JAVAHOME,javahome, classpath, LDLIBRARY,ldpath)
