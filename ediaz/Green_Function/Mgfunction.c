/* Green Function modelling
   from analytical formula for Isotropic homogeneus media.
 
  u_{i}(\bar{x},t) = \frac{1}{4\pi\rho}(3\gamma_i \gamma_j)
     \int_{R/\alpha}ˆ{R/\beta} \tau X_{0}(t-\tau)d\tau 
     + \frac{1}{4\pi\rho\alpha^2} \gamma_i \gamma_j \frac{1}{R} X_0 (t-R/\alpha)
     + \frac{1}{4\pi\rho\beta^2}(\delta_{ij} -\gamma_i \gamma_j)\frac{1}{R} X_0 (t - \frac{R}{\beta})

 Units of R, alpha, beta, rho must be consistant. Is up to you (the user) to get the seismogram
 rigth
 
 R     : (float) [100]  distance from source to receiver
 theta : (float) [0]    angle between receiver line and x1 axis
 alpha : (float) [2000] P-wave velocity
 Beta  : (float) [1000] S-wave velocity
 rho   : (float) [2500] density
 force : (int)   [1]    force direction (1  or 2 ) 
 stdin : (file of floats) source function 
 
  -----------------> x1
 |\ ) \theta
 | \ 
 |  \
 |   \
 |    \     medium properties: alpha, beta, rho
 |	   R
 |
 |
 |	 
 x3
 
 Written by Esteban D\'{i}az, Fall 2011 
 for Introduction to Seismology

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
    //time sampling to get from source file
    int nt,inttau;
	float ot,dt,t,tau,x;
	int i1;	
	float snf,spff,ssff,pi;
	
	float nf,ffp,ffs;
	float lo,hi;
	
	float r,theta; //receiver variables
	float alpha, beta , rho; //properties of the medium  
	float gamma1,gamma2;
    int force; //force direction :1 parallel to x1 axis, 2 par to x2.
	// I/O files 2D wfield (out), source function (in) 
    float **trace, *source;
	sf_file  in, out;

	//================= End of variable declaration =====================
	
    sf_init (argc,argv);
    out = sf_output("out");
	in = sf_input ("in"); 
	
	//====  Read parameters from command line ==================	

	 
	//R: source-receiver distance?
	if(! sf_getfloat("R",&r)) r=100.0; 
    //theta: angle between receiver line and x1 axis? 
	if(! sf_getfloat("theta",&theta)) theta=0.0; 
	//alpha: P-wave velocity
	if(! sf_getfloat("alpha",&alpha)) alpha=2000.0;
	//beta: S-wave velocity
	if(! sf_getfloat("beta",&beta)) beta=1000.0;
	//rho : density
	if(! sf_getfloat("rho",&rho)) rho=2500.0;	
	//force: 
	if(! sf_getint("force",&force)) force=1;	
 

	//======== Get parameters from input file  =================
	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input");
	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d3= in input");
	
    trace = sf_floatalloc2(2,nt);
    source= sf_floatalloc(nt);

	
	sf_floatread(source,nt,in);
    sf_floatwrite(source,nt,out);

    pi=acos(-1);

	snf=1/(4*pi*rho);
	spff=snf/(alpha*alpha*r);
	ssff=snf/(beta*beta*r);
	
	snf*=1/(pow(r,3));
	
	lo=r/alpha;
	hi=r/beta;
    
	gamma1= cos(theta*pi/180);
	gamma2= sin(theta*pi/180);
	
	
	for (i1=0; i1<nt; i1++) {
		nf=0;ffp=0;ffs=0;
		t=ot+dt*i1;
		//Near field:
		nf=0;
		tau=lo;
		while (tau<hi+0.001) {
			inttau= floorf((t-tau)/dt);
			x=0;
			if (inttau>=0) x=source[inttau];			
			nf += tau*x*dt ;
			tau +=dt;
		}
		nf*=snf;
		
		inttau=floorf((t-lo)/dt);
		if(inttau >=0)  ffp=source[inttau];
        ffp*=spff;
		
		inttau=floorf((t-hi)/dt);
		if(inttau >=0)  ffp=source[inttau];
		ffs*=ssff;

		trace[1,i1]=nf+ffp+ffs;
		trace[2,i1]=nf+ffp+ffs;

		sf_floatwrite(trace[1,i1],1,out);
		sf_floatwrite(trace[2,i1],1,out);		
		
		
	}
	
	
	exit(0);
}