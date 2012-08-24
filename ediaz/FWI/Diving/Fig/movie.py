#!/usr/bin/env python

import commands as cm



list_vpl=cm.getoutput("ls -1 sgrad_s-*.pdf |sed -e \"s/sgrad//g\"|sed -e \"s/.pdf//g\" |awk '{printf(\"%s \",$1)}' ")


vpls=list_vpl.split()

for filename in vpls:

    var=str(filename)

    a='%--------'+'''
\\begin{frame} 
    \\begin{columns}
        \\column{0.5\\textwidth}
            \\plot{sgrad%s}{width=\\textwidth}{\\klabellarge{10}{50}{source gradient}}
            \\plot{pen-cig%s}{width=0.85\\textwidth}{\\klabellarge{10}{-10}{\\blue{time-lag}}\\klabellarge{8}{-20}{\\red{penalized}}}
        \\column{0.5\\textwidth}
            \\plot{rgrad%s}{width=\\textwidth}{\\klabellarge{10}{50}{receiver gradient}}
            \\plot{adj_srcs%s}{width=0.85\\textwidth}{\\klabellarge{10}{-10}{\\blue{source adjoint src}}\\klabellarge{8}{-20}{\\red{rec adjoint src}}}
    \\end{columns}
\\end{frame} 
'''%(var,var,var,var)+'%--------'

    print a 
