/* Second derivative along third axis

 */
/*
  Copyright (C) 2011 Colorado School of Mines

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>

int main (int argc, char* argv[])
{
    int n1,n2,i2,i1,n3, n;
    float *uo,*u1,*u2,*der;
    float c0,c1,c2,dt;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d3",&dt)) sf_error("No d3= in input");

    n3 = sf_leftsize(in,2);

    uo  = sf_floatalloc(n1*n2);
    u1  = sf_floatalloc(n1*n2);
    u2  = sf_floatalloc(n1*n2);
    der = sf_floatalloc(n1*n2);


    c0=(1/dt)*(1/dt);
    c1=-2*c0;
    c2=c0;

    for (i1=0; i1<n1*n2; i1++){
	uo[i1]=0.0;
	u1[i1]=0.0;
	u2[i1]=0.0; 
	der[i1]=0.0;	
    }



    sf_floatread(u1,n1*n2,in);
    sf_floatread(u2,n1*n2,in);
    
    for (i1=0; i1<n1*n2; i1++){
	der[i1]=c1*u1[i1]+c2*u2[i1];	
    }
    sf_floatwrite(der,n1*n2,out);



    for (i2=1; i2 < n3-1; i2++) {

	for (i1=0; i1<n1*n2; i1++){
	   uo[i1]= u1[i1];
           u1[i1]= u2[i1];
        } 	
        sf_floatread(u2,n1*n2,in);
        for (i1=0; i1<n1*n2; i1++){ 
     	   der[i1]=c0*uo[i1] +c1*u1[i1]+c2*u2[i1];	
	}		
        sf_floatwrite(der,n1*n2,out);
    }

    for (i1=0; i1<n1*n2; i1++){
	der[i1]=c1*u2[i1]+c0*u1[i1];	
    }
    sf_floatwrite(der,n1*n2,out);

    exit(0);
}
