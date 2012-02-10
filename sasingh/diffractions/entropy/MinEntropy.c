/* Measure of diffraction focusing using Negentropy */
/*
 	Copyright (C) 2008 Colorado School of Mines
 	
 	This program is free software; you can redistribute it and/or modify
 	it under the terms of the GNU General Public License as published by
 	the Free Software Foundation; either version 2 of the License, or
 	(at your option) any later version.
 	
 	This program is distributed in the hope that it will be useful,
 	but WITHOUT ANY WARRANTY; without even the implied warranty of
 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 	GNU General Public License for more details.
 	
 	You should have received a copy of the GNU General Public License
 	along with this program; if not, write to the Free Software
 	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 	*/

#include <rsf.h>
#include <math.h>
#include <stdio.h>
/*
#ifdef _OPENMP
#include <omp.h>
#endif
*/

int main(int argc, char* argv[])
{
//    int ompnth=1;

    sf_file in, out; /* Input and output files */
    sf_axis az,ax,a3;
    int nz,nx,n3; 
    int box1,box2,klo1,khi1,klo2,khi2,kmid2,kmid1;
    int ix,iz,i3;

    float *sumG,***d,***dat,***neg,***dN;
    float neg1,AmpNorm,h,h1; 
    double scalea, logs;
    
/*---------------------------------------------------------*/
    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in  = sf_input("in");
    /* standard output */
    out = sf_output("out");
/*
#ifdef _OPENMP
   ompnth = omp_init();
#endif
*/    
    /* parameters from input file*/
    az=sf_iaxa(in,1); sf_setlabel(az,"z"); nz = sf_n(az);
    ax=sf_iaxa(in,2); sf_setlabel(ax,"x"); nx = sf_n(ax);
 	a3=sf_iaxa(in,3); sf_setlabel(a3,"y"); n3 = sf_n(a3);


    /* parameter from the command line (i.e. box1=50 box2=50 ) */
    if (!sf_getint("box1",&box1)) sf_error("Need box1=");
    if (!sf_getint("box2",&box2)) sf_error("Need box2=");
    
    /* allocate floating point array */
    dat = sf_floatalloc3 (nz,nx,n3);
    
    /* initialise the size of the searching box*/
    int s1= nz-(box1);
    int s2= nx-(box2);
    int bm1=box1/2;
    int bm2=box2/2;

    /*initialise the mid-point of each box) */

    sumG =sf_floatalloc  (n3);
    d    =sf_floatalloc3 (nz,nx,n3);
    dN   =sf_floatalloc3 (nz,nx,n3);
    neg  =sf_floatalloc3 (nz,nx,n3);

    sf_floatread(dat[0][0],nz*nx*n3,in); 

// Global Sum
    for     (i3=0 ; i3<n3; ++i3){
            sumG[i3]=0;
      for   (ix=0; ix<nx; ++ix){
        for (iz=0; iz<nz; ++iz){
	         d[i3][ix][iz] = dat[i3][ix][iz] * dat[i3][ix][iz]; //make all amplitudes positive
           sumG[i3] += d[i3][ix][iz];
        }
	  }
    }

/*
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(klo2,klo1,kmid2,kmid1,khi1,khi2,neg1,AmpNorm,scalea,logs,h1,h) 
#endif
*/


//initialise the box for negentropy
  for     (i3=0  ; i3<n3    ; ++i3  ){ 
          klo2=0;klo1=0;

    for   (klo2=0; klo2<(s2); ++klo2){
      for (klo1=0; klo1<(s1); ++klo1){
	
//intialise parameters
          neg1=0.0;scalea=0.0;logs=0.0;
          AmpNorm=0.0;h=0,h1=0;

//Set upper limit of searching box and midpoints   
          khi2=klo2+box2; khi1=klo1+box1;
          kmid2=klo2+bm2; kmid1=klo1+bm1;        

//Sum values in each box         
         for    (ix=klo2; ix<khi2; ++ix){
            for  (iz=klo1; iz<khi1; ++iz){
                 AmpNorm       =(d[i3][ix][iz])/(sumG[i3]);
                 dN[i3][ix][iz]=AmpNorm;

//Gaussian operator
                 h =(((iz  -kmid1) * (iz  -kmid1)) +       \
                     ((ix  -kmid2) * (ix  -kmid2))) / (2*box1*box2*1.41);
                 h1=(((box1*0.5  ) * (box1*0.5  )) +       \
                     ((box2*0.5  ) * (box2*0.5  )))/ (2*box1*box2*1.41);
                 h =exp(-4*h );
                 h1=exp(-4*h1);

                 if (h1>= h) scalea=0;
                    else
                      scalea =   AmpNorm * (h-h1);
                 if (AmpNorm==0) logs= 0;
		            else 	     logs= scalea*scalea;

/* logs can be different functions:
            		else {logs=log(scalea);}
		            logs=log(scalea); */

			        neg1 += (scalea*logs);
	        }
          }
          neg[i3][kmid2][kmid1] = neg1/(box1*box2);
        }
      }
    }
    sf_floatwrite(neg[0][0]  ,nz*nx*n3 ,out);  //write negentropy
  

    exit(0);
}
