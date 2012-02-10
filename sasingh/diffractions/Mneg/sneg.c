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

int main(int argc, char* argv[])
{
 
    int n1, n2, box1,box2,klo1,khi1,klo2,khi2,kmid1,kmid2;
    int o1,ix,iz,m1;

    float  **d,**dat,**neg,**dN;
    float sumG,neg1,AmpNorm,h,h1; 
    double scalea, logs;
    sf_file in, out; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in  = sf_input("in");
    /* standard output */
    out = sf_output("out");
    
    /* parameters from input file*/
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

    /* parameter from the command line (i.e. box1=50 box2=50 ) */
    if (!sf_getint("box1",&box1)) sf_error("Need box1=");
    if (!sf_getint("box2",&box2)) sf_error("Need box2=");
    
    /* allocate floating point array */
    dat = sf_floatalloc2 (n1,n2);
    
    /* initialise the size of the searching box*/
    klo1=0;
    khi1=box1;
    klo2=0;
    khi2=box2;
    int s1= n1-(box1);
    int s2= n2-(box2);

    /*initialise the mid-point of each box) */
    kmid1=((klo1+khi1)*0.5);
    kmid2=((klo2+khi2)*0.5)-1;

    d    =sf_floatalloc2 (n1,n2);
    dN   =sf_floatalloc2 (n1,n2);
    neg  =sf_floatalloc2 (n1,n2);
    for (ix=0; ix<n2; ix++) {
      for (iz=0; iz<n1; iz++) {
        neg[ix][iz]=0.0;
        dat[ix][iz]=0.0;
        d[ix][iz]=0.0;
	      dN[ix][iz]=0.0;
      }
    }
    sf_floatread(dat[0],n1*n2,in); 
    sumG=0;

// Global Sum
    for (ix=o1; ix<n2; ++ix){
      for (iz=o1; iz<n1; ++iz){
	        d[ix][iz]= dat[ix][iz] * dat[ix][iz]; //make all amplitudes positive
          sumG += d[ix][iz];// sum of each trace
      }
	  }
   


//initialise the box for negentropy
   for (klo2=0; klo2<(s2) ; ++klo2) {
       kmid2++; 
       m1=kmid1-1;
      for (klo1=0; klo1<(s1); ++klo1){
         m1++;
	
//intialise parameters
    neg1=0.0;scalea=0.0;logs=0.0;AmpNorm=0.0;h=0,h1=0;

//Set upper limit of searching box    
    khi2=klo2+box2; khi1=klo1+box1;
         
//Sum values in each box         
         for (ix=klo2; ix<khi2; ++ix){
            for (iz=klo1; iz<khi1; ++iz){
                 AmpNorm =(d[ix][iz])/(sumG);
                 dN[ix][iz]=AmpNorm;

//Penalty operator
                 h=(((iz-m1)*(iz-m1))+((ix-kmid2)*(ix-kmid2)))/(2*box1*box2*1.41);
                 h=exp(-4*h);

                 h1=(((box1*0.5)*(box1*0.5))+((box2*0.5)*(box2*0.5)))/(2*box1*box2*1.41);
                 h1=exp(-4*h1);

                 if (h1>= h) scalea=0;
                    else
                      scalea=box1*box2*AmpNorm*(h-h1);
                 if (AmpNorm==0) {logs=0;}
		                else {	logs= scalea*scalea;}

/* logs can be different functions:
            		else {logs=log(scalea);}
		            logs=log(scalea); */

			           neg1 += (scalea*logs);
	        	}
         }
        	       neg[kmid2][m1] = neg1/(box1*box2);
      }
    }

         sf_floatwrite(neg[0]  ,n1*n2 ,out);  //write negentropy

    exit(0);
}
