#!/usr/bin/env python
"""
Nonlinear Gradient-based solver
"""

import os, sys, shutil, glob, commands, string, time, math, datetime
import rsf.user.ivlad as ivlad
from string import lower

try:
    import rsf
except:
    import rsfbak as rsf

class Par:

    def __init__(self,argv=sys.argv):
        self.noArrays = True
        self.prog = argv[0]
        self.__args = self.__argvlist2dict(argv[1:])

    def __argvlist2dict(self,argv):
        """Eliminates duplicates in argv and makes it a dictionary"""
        argv = self.__filter_equal_sign(argv)
        args = {}
        for a in argv:
            key = a.split('=')[0]
            args[key] = a.replace(key+'=','')
        return args

    def __filter_equal_sign(self,argv):
        """Eliminates "par = val", "par= val" and "par =val" mistakes."""
        argv2 = []
        # Could not use the simpler 'for elem in argv'...argv.remove because
        # lonely '=' signs are treated weirdly. Many things did not work as
        # expected -- hence long and ugly code. Test everything.
        for i in range( len(argv) ):
            if argv[i] != '=':
                if argv[i].find('=') != 0:
                    if argv[i].find('=') != -1:
                        if argv[i].find('=') != len(argv[i])-1:
                            argv2.append(argv[i])
        return argv2

    def __get(self, key, default):
        """Obtains value of argument from dictionary"""
        if self.__args.has_key(key):
            return self.__args[key]
        elif str(default):
            return default
        else:
            sys.stderr.write( '\n  Argument %s= not given to %s \n' %
                              (key, self.prog))
            sys.exit(1)

    def string(self, key):
        """Returns string argument given to program"""
        return self.__get(key, default=None)

    def int(self,key,default=None):
        """Returns integer argument given to program"""
        try:
            return int( self.__get(key, default) )
        except:
            sys.stderr.write( '\n  Argument %s= to %s must be integer\n' %
                              (key, self.prog))
            sys.exit(1)

    def float(self,key,default=None):
        """Returns float argument given to program"""
        try:
            return float( self.__get(key, default) )
        except:
            sys.stderr.write( '\n  Argument %s= to %s must be float\n' %
                              (key, self.prog))
            sys.exit(1)

    def bool(self,key,default=None):
        """Returns bool argument given to program"""
        val = self.__get(key, default)
        val = lower(str(val)) # No, combining with line above does not work
        if val == 'y' or val == 'true':
            return True
        elif val =='n' or val == 'false':
            return False
        else:
            msg = ('\n  Argument %s= to %s must be bool (y/n, True/False) \n' %
                   (key, self.prog))
            sys.stderr.write(msg)
            sys.exit(1)


####################################################################

def write2log(logfile, statement):
    if isinstance(logfile,file):
        logfile.write(statement)
    elif isinstance(logfile,str):
        handle = open(logfile,'a')
        handle.write(statement+'\n')
        handle.close()

####################################################################

def log_init(log_file):
    global logfile
    logfile = log_file
    handle = open(logfile,'a').close() # Erase log

####################################################################

def report_init(report_file):
    global repfile
    repfile = report_file
    handle = open(repfile,'a').close() # Erase log

####################################################################

def warning(message):
    'Screen output and logging message '

    global logfile
    print message 

    write2log(logfile, message)

####################################################################

def report_warning(message):
    'Logging message '

    global repfile

    write2log(repfile, message)

#############################################################################

def clean_make_outdir(outdir):
    'Clean outdir with sfrm to avoid disk leaks'

    if os.path.exists(outdir):
        # clean up invalid files which will make sfrm fail
        invalids = commands.getoutput('sfinvalid dir=' + outdir)
        if invalids.strip() != '':
            invalist = invalids.split()
            for ifile in invalist:
                os.remove(ifile)
        rsf_files =  glob.glob(os.path.join(outdir, '*' + ivlad.ext))
        for rsf_file in rsf_files:
            os.system('sfrm ' + rsf_file)
        shutil.rmtree(outdir)
    os.mkdir(outdir)

#############################################################################

def clean_outdir(outdir):
    'Clean outdir with sfrm to avoid disk leaks'

    if os.path.exists(outdir):
        # clean up invalid files which will make sfrm fail
        invalids = commands.getoutput('sfinvalid dir=' + outdir)
        if invalids.strip() != '':
            invalist = invalids.split()
            for ifile in invalist:
                os.remove(ifile)
        rsf_files =  glob.glob(os.path.join(outdir, '*' + ivlad.ext))
        for rsf_file in rsf_files:
            os.system('sfrm ' + rsf_file)
        shutil.rmtree(outdir)

#############################################################################

def get_concatenation_axis(filenm):
    'Calculates the last dim of filenm, adds one, converts to string'

    cmd = '<%s sffiledims parform=n' % filenm
    nds = int(commands.getoutput(cmd).split(':')[0])
    return str(nds+1)

####################################################################

def execute(command):
    'Echo command to screen, write it to log, execute it'
    global logfile
    print command + '\n' 

    write2log(logfile, command+'\n')
    os.system(command)

####################################################################

def add_zeros(i, n):
    'Generates string of zeros + str(i)'

    ndigits_n = int(math.floor(math.log10(n-1)))
    if i == 0:
        nzeros = ndigits_n
    else:
        nzeros = ndigits_n - int(math.floor(math.log10(i)))

    return nzeros*'0'+str(i)

####################################################################

def getnorm(file):
    'Squares of L2 norm of input'

    norm = 0
    norm = float(commands.getoutput('<'+file+' sfnorm '))

    return norm

#############################################################################

def getinnerprod(file1,file2):
    'Get innerproduct between two files'

    result = 0
    result += float(commands.getoutput('sfinnerprod <' + file1 + ' y=' + file2))

    return result

#############################################################################

def model_update(x_old, a, delta_x, x_new):
    'Update the model with given direction and step length'

    # x_new = x_old + a*delta_x
    execute('<%s sfscale axis=123 | sfaxplusy a=%1.15E y=%s >%s' % (delta_x.lstrip(),a,x_old,x_new))

#############################################################################

def conj_update(x_old, a, delta_x, x_new):
    'Update the conjugate gradient and normalize'

    # x_new = x_old + a*delta_x   
    execute('<%s sfaxplusy a=%1.15E y=%s >%s' % (delta_x.lstrip(),a,x_old,x_new))

#############################################################################

def optimized_alpha(objt, alpha):
    ' Find the vertex of the objective function curve (hyperbola approximation) '

    A = 0
    B = 0
    tmp = ((alpha[2]*alpha[2])/(alpha[1]*alpha[1]))

    B = (objt[2] - objt[0] - tmp*(objt[1]-objt[0]))/(alpha[2] - tmp*alpha[1])
    A = (objt[1] - objt[0] - alpha[1]*B)/(alpha[1]*alpha[1])
    vertex = -0.5*B/A

    return vertex

#############################################################################
 
def str2bool(v):
  return v.lower() in ("y", "true", "yes", "1")

#############################################################################

def main(argv=sys.argv):

    if not os.path.isfile('runsolver.par'):
        print 'sover parameters file is missing! Please check it!'
        return unix_error
    else:
        par_file = open('runsolver.par','r')
        par = par_file.readlines()

    objt = par[0].strip()
    objt_args = par[1].strip()
    objt_uses_tag = str2bool(par[2].strip())

    grad = par[3].strip()
    grad_args = par[4].strip()
    grad_uses_tag = str2bool(par[5].strip())
    
    mod = par[6].strip()
    sol = par[7].strip()
    movie = par[8].strip()

    niter = string.atoi(par[9])
    nstep = string.atoi(par[10])
    dstep = string.atof(par[11])

    cg_tag = str2bool(par[12].strip())

    if os.path.isfile(sol):
        print 'Target %s' %sol + ' has been build, nothing to do!'
        return ivlad.unix_success

    print 'read parameters done, please double check:\n' 
    print 'objt=',objt
    print 'objt_args=',objt_args
    print 'objt_uses_tag=',objt_uses_tag
    print 'grad=',grad
    print 'grad_args=',grad_args
    print 'grad_uses_tag=',grad_uses_tag
    print 'mod=',mod
    print 'sol=',sol
    print 'movie=',movie
    print 'niter=',niter
    print 'nstep=',nstep
    print 'dstep=',dstep
    print 'cg=',cg_tag
    print '\n'

    tolerance = 1e-20             # criterion for convergence  
    objt_inp = []
    objt_out = []
    grad_inp = []
    grad_out = []

    # create input output symbols for operator computing objective function and gradient
    if objt_uses_tag:
        objt_inp=(lambda x: ' --input='+x)
        objt_out=(lambda x: ' --output='+x)
    else:
        objt_inp=(lambda x: ' <'+x)
        objt_out=(lambda x: ' >'+x)

    if grad_uses_tag:
        grad_inp=(lambda x: ' --input='+x)
        grad_out=(lambda x: ' --output='+x)
    else:
        grad_inp=(lambda x: ' <'+x)
        grad_out=(lambda x: ' >'+x)

    objt_file = []
    grad_file = []
    conj_file = []
    x_file    = []

    sol_nameroot = os.path.splitext(os.path.basename(sol))[0]

    # initialize log file
    log = 'log_'+sol_nameroot+'.asc'
    log_init(log)

    # initialize inversion log file
    report = 'report_'+sol_nameroot+'.asc'
    report_init(report)

    ###############   pre-processing, determine what action should take   #############


    status_dir = './status_' + sol_nameroot + '/'

    if not os.path.isdir(status_dir):
        warning('\n')
        warning('**********     Begin Nonlinear Gradient-based solver at: ' + time.ctime() + '    **********\n')
        print '**********     No previous inversion exits!     **********\n'
        os.mkdir(status_dir)

        iter_old = 0
        inversion_flag = 0
        line_flag = 0
        
        iter_log = open(status_dir+'iter_log.asc','w')
        iter_log.write(str(iter_old)+'\n')
        iter_log.write(str(inversion_flag)+'\n')
        iter_log.write(str(line_flag)+'\n')
        iter_log.close()   

        iter_string = 'iter%03d' %0 

        command = 'sfcp ' + mod + ' ' + status_dir + 'sol'+'-'+iter_string+'.rsf'
        os.system(command)

    else:
        print '**********     Previous inversion exits!     **********\n'
 
        iter_log = open(status_dir+'iter_log.asc','r')
        info = iter_log.readlines()
        iter_old = string.atoi(info[0])
        inversion_flag = string.atoi(info[1])
        line_flag = string.atoi(info[2])
 
        iter_new  = iter_old + 1
        iter_last = iter_old - 1

        iter_string      = 'iter%03d' % iter_old
        iter_string_new  = 'iter%03d' % iter_new
        iter_string_last = 'iter%03d' % iter_last

        ls_string   = '-T%02d' %line_flag

        if inversion_flag == 4:
            print 'Target %s' %sol + ' has been build, nothing to do!'
            return unix_error

        if inversion_flag == 3:
            if os.path.isfile(status_dir+'sol'+'-'+iter_string_new+'.rsf'):
              print '**********     Previous computation: model update   **********\n'

              if iter_new == niter:
                inversion_flag = 4
              else :
                inversion_flag = 0
                iter_old = iter_new
                iter_string = 'iter%03d' % iter_old
 
                line_flag = 0

        if inversion_flag == 2:

            if line_flag == 2:
              if os.path.isfile(status_dir+'objt'+'-'+iter_string+ls_string+'.rsf'):
                print '**********     Previous computation: line search #%d   **********\n' %line_flag
                inversion_flag = 3

            if line_flag == 1:
              if os.path.isfile(status_dir+'objt'+'-'+iter_string+ls_string+'.rsf'):
                print '**********     Previous computation: line search #%d   **********\n' %line_flag
                line_flag += 1

        if inversion_flag == 1:
            if os.path.isfile(status_dir+'grad'+'-'+iter_string+'.rsf'):
              print '**********     Previous computation: gradient   **********\n'
              inversion_flag += 1
              line_flag = 1

        if inversion_flag == 0:
            objt_file.append(status_dir+'objt'+'-'+iter_string+'.rsf')

            if os.path.isfile(objt_file[0]):
              print '**********     Previous computation: objective function   **********\n'

              # record the value of f(x_i) in file
              frame2 = status_dir + 'objt_movie'+'-'+iter_string+'.rsf'
              command2 = ' sfmath n1=1 o1=0 d1=1 output="%g" > ' %getnorm(objt_file[0]) + frame2
              os.system(command2)

              if getnorm(objt_file[0]) < tolerance: 
                inversion_flag = 4
              else:
                inversion_flag += 1

        info[0] = str(iter_old)+'\n'
        info[1] = str(inversion_flag)+'\n'
        info[2] = str(line_flag)+'\n'

        iter_log = open(status_dir+'iter_log.asc','w')
        iter_log.writelines(info)
        iter_log.close()
    #############################   pre-processing done   ############################


    ########   processing, take action depending on pre-processing parameter   ########
    
    if inversion_flag == 0: # compute O.F.
        warning('********************************************************************************************************' )
        warning('**********     compute objective function for iteration %d of %d at: ' %(iter_old,niter) + time.ctime() + ' **********')
        warning('********************************************************************************************************\n' )

        # compute objective function f(x_i)
        x_file.append(status_dir+'sol'+'-'+iter_string+'.rsf')
        objt_file.append(status_dir+'objt'+'-'+iter_string+'.rsf')

        command = objt + ' ' + objt_args + ' ' + objt_inp(x_file[0]) + ' --itag=' + iter_string + objt_out(objt_file[0]) 
        execute(command)

        return ivlad.unix_success

    if inversion_flag == 1: # compute gradient
        warning('***********************************************************************************************' )
        warning('**********     compute gradient for iteration %d of %d at: '%(iter_old,niter) + time.ctime() + '  **********')
        warning('***********************************************************************************************\n' )

        # compute grad f(x_i)
        grad_file.append(status_dir+'grad'+'-'+iter_string+'.rsf')
        objt_file.append(status_dir+'objt'+'-'+iter_string+'.rsf')
        x_file.append(status_dir+'sol'+'-'+iter_string+'.rsf')

        command = grad + grad_inp(objt_file[0]) + ' ' + grad_args + ' --itag=' + iter_string + ' --model=' + x_file[0] + ' ' + grad_out(grad_file[0])
        execute(command)

        return ivlad.unix_success

    if inversion_flag == 2: # compute O.F. in line search
        warning('****************************************************************************************************' )
        warning('**********     compute line search #%d for iteration %d of %d at: '%(line_flag,iter_old,niter) + time.ctime() + '  *********')
        warning('****************************************************************************************************\n' )

        ls_string = '-T%02d' % line_flag

        # compute objective function f(x_i)
        x_file.append(status_dir+'sol'+'-'+iter_string+'.rsf')
        x_file.append(status_dir+'sol'+'-'+iter_string+ls_string+'.rsf')

        objt_file.append(status_dir+'objt'+'-'+iter_string+ls_string+'.rsf')
        grad_file.append(status_dir+'grad'+'-'+iter_string+'.rsf')

        if cg_tag:
            print 'Use conjugate gradient direction!'
            if iter_old == 0:
                conj_file.append(status_dir+'conj'+'-'+iter_string+'.rsf')
                command = '< ' + grad_file[0] + ' sfwindow squeeze=n > ' + conj_file[0]
                os.system(command)

                model_update(x_file[0],dstep*line_flag,conj_file[0],x_file[1])

            else:
                # notice that here [1] means the file in the previous iteration, while [0] means the file
                # in the current iteration
                grad_file.append(status_dir+'grad'+'-'+iter_string_last+'.rsf')

                conj_file.append(status_dir+'conj'+'-'+iter_string+'.rsf')
                conj_file.append(status_dir+'conj'+'-'+iter_string_last+'.rsf')

                beta = getnorm(grad_file[0])/getnorm(grad_file[1])
                conj_update(grad_file[0],beta,conj_file[1],conj_file[0])

                model_update(x_file[0],dstep*line_flag,conj_file[0],x_file[1])
        else:
            print 'Use steepest descent direction!'
            model_update(x_file[0],dstep*line_flag,grad_file[0],x_file[1])

        command = objt + objt_inp(x_file[1]) + ' ' + objt_args + ' --itag=' + iter_string + ls_string + objt_out(objt_file[0])       
        execute(command)

        return ivlad.unix_success

    if inversion_flag == 3: # compute alpha.
        warning('***********************************************************************' )
        warning('**********     compute step length at iteration %d of %d   **********' %(iter_old,niter))
        warning('***********************************************************************\n' )

        x_file.append(status_dir+'sol'+'-'+iter_string+'.rsf')
        x_file.append(status_dir+'sol'+'-'+iter_string_new+'.rsf')
        grad_file.append(status_dir+'grad'+'-'+iter_string+'.rsf')
        conj_file.append(status_dir+'conj'+'-'+iter_string+'.rsf')

        # linear search for step length 
        objt_trial = [0]*nstep 
        alpha = [0]*nstep

        for j in range(nstep):
           ls_string = '-T%02d' %j
           alpha[j] = j * dstep
           if j == 0:
             objt_trial[j] = getnorm(status_dir+'objt'+'-'+iter_string+'.rsf')
           else:
             objt_trial[j] = getnorm(status_dir+'objt'+'-'+iter_string+ls_string+'.rsf')

        # update model
        best_alpha = optimized_alpha(objt_trial,alpha)
        #if best_alpha < 0:
            #best_alpha = alpha[objt_trial.index(min(objt_trial))]
        #    best_alpha = 0.5*dstep
        if best_alpha > nstep*dstep:
             best_alpha = alpha[objt_trial.index(min(objt_trial))]
             if best_alpha == 0:
               best_alpha = 0.5*dstep

        warning('the value for optimized alpha is %g\n' %(best_alpha))
        report_warning('iteration %d: (%g,%g), (%g,%g), (%g,%g) --> computed alpha: %g, used alpha: %g\n' %(iter_old,alpha[0],objt_trial[0],alpha[1],objt_trial[1],alpha[2],objt_trial[2],optimized_alpha(objt_trial,alpha),best_alpha))

        if cg_tag:
            model_update(x_file[0],best_alpha,conj_file[0],x_file[1])
        else:
            model_update(x_file[0],best_alpha,grad_file[0],x_file[1])

        # clear all the temporary files and directories
        os.system('sfrm ' + iter_string + '*.rsf')

    if inversion_flag == 4: # clean up

        iter_zero_string = 'iter%03d' %0 

        x_file.append(status_dir+'sol'+'-'+iter_string+'.rsf')

        # Return solution
        os.system('sfcp ' + x_file[0] + ' ' + sol)

        cdim_movie  = get_concatenation_axis(mod)
        cdim_movie1 = str(4)
        cdim_movie2 = str(2)
        cdim_movie3 = get_concatenation_axis(status_dir+'objt-'+iter_zero_string+'.rsf')
        command  = 'sfcat axis=' + cdim_movie  + ' '
        command1 = 'sfcat axis=' + cdim_movie1 + ' '
        command2 = 'sfcat axis=' + cdim_movie2 + ' '
        command3 = 'sfcat axis=' + cdim_movie3 + ' '

        for i in range(niter):
            iter_string = '-iter%03d' %i

            command  += status_dir+'sol'+iter_string+'.rsf '
            command1 += status_dir+'grad'+iter_string+'.rsf '
            command2 += status_dir+'objt_movie'+iter_string+'.rsf '
            command3 += status_dir+'objt'+iter_string+'.rsf '

        movie_grad = 'movie_grad_'+sol_nameroot+'.rsf'
        movie_resi = 'movie_resi_'+sol_nameroot+'.rsf'
        movie_objt = 'curve_objt_'+sol_nameroot+'.rsf'

        os.system(command  + ' > ' + movie )
        os.system(command1 + ' | sfwindow > ' + movie_grad )
        os.system(command2 + ' | sfwindow > ' + movie_objt )
        os.system(command3 + ' > ' + movie_resi )

        # clear all the temporary files and directories
        for i in range(niter):
            iter_string = '-iter%03d' %i
            os.system('rm pbs'+iter_string+'* -rf')

        warning('*****************************************************\n' )
        warning('**********     Solver running details are recorded in: %s     **********\n' %log)
        warning('**********     The solution file is: %s     **********\n' %sol)
        warning('**********     The movie for solution is: %s     **********\n' %movie)
        warning('**********     The movie for gradient is: %s     **********\n' %movie_grad)
        warning('**********     The movie for residual is: %s     **********\n' %movie_resi)
        warning('**********     The objective function curve is: %s     **********\n' %movie_objt)
        warning('**********     inversion finished!   **********\n')
        warning('**********     end time: ' + time.ctime() + '     **********\n')
        warning('*****************************************************\n' )

        return ivlad.unix_success

    os.system('./Nonlinear_CG_solver.py')

##############################################

if __name__ == '__main__':

    try:
        status = main()
    except:
        status = ivlad.unix_error

    sys.exit(status)

