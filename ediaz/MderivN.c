/* First derivative with a maximally linear FIR differentiator. 

Division by dx is not included. */
/*
  Copyright (C) 2004 University of Texas at Austin

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
    int n1,n2, i1,i2,i3, n,n3,ndim;
    float c0, c1, c2;
    float dt,d3;
    float **um,**uo,**ut, **der;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");     /* Wavefield in*/
    out = sf_output("out");  /* Wavefield 2nd derivative out */



   /* Wavefield second derivative on third axis */

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d3",&d3)) sf_error("No d3= in input");

    um=sf_floatalloc2(n1,n2);
    uo=sf_floatalloc2(n1,n2);
    ut=sf_floatalloc2(n1,n2);
    der=sf_floatalloc2(n1,n2);


    for    (i1=0; i1<n1; i1++) {
	for(i2=0; i2<n2; i2++) {
	    um[i1][i2]=0;
	    uo[i1][i2]=0;
	    ut[i1][i2]=0;    
           der[i1][i2]=0;
	}
    }


    c0 = 1.0; c1=-2.0; c2=1.0;
  
    c0= c0/(dt*dt); c1=c1/(dt*dt); c2=c2/(dt*dt);
   
    /* Filter order */

 
    sf_floatread(uo[0],n1*n2,in);
    sf_floatread(ut[0],n1*n2,in);
    fprintf(stderr,"%4d%4d\n",n1,n2);
    for    (i1=0; i1<n1; i1++) {
        for(i2=0; i2<n2; i2++) {
	   der[i1][i2]= -2*uo[i1][i2]+ut[i1][i2];
           if (uo[i1][i2] >0.0) fprintf(stderr, "hola\n");
           fprintf(stderr,"%f%f%f\n",der[i1][i2],uo[i1][i2],ut[i1][i2]);
        }
    }    
    um=uo; uo=ut;
    sf_floatwrite(der[0],n1*n2,out);

   /*

    for (i3=1; i3 < n3 -1 ; i3++) {
        sf_floatread(ut[0],n1*n2,in);
        for    (i1=0; i1<n1; i1++) {
	    for(i2=0; i2<n2; i2++) {               
	      der[i1][i2]= um[i1][i2]-2*uo[i1][i2]+ut[i1][i2];
	    }
        }    
        um=uo; uo=ut;
        sf_floatwrite(der[0],n1*n2,out);
    }

    for    (i1=0; i1<n1; i1++) {
        for(i2=0; i2<n2; i2++) {
	   der[i1][i2]= um[i1][i2]-2*uo[i1][i2];
        }
    }    
    sf_floatwrite(der[0],n1*n2,out);


    /*Deallocate **um,**uo,**ut, **der; */

    exit(0);
}
