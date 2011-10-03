/*Acoustic time-domain FD modeling, automatically determines to use 2D or 3D modeling based on velocity file structure.
 *
 * VERSION 1.0 - May 6, 2011

Acoustic wave equation finite difference modeling in both 2D and 3D using
an explicit time-domain solver.

**** PLEASE SEE THE SCONSTRUCT in book/tutorial/awe for examples on how to 
use this program ****

This program is designed to be as generic as possible, and allows you to use files
with arbitrary models, arbitrary source and receiver geometries, and arbitrary
source types (acceleration or displacement).

========= OPTIONS ============
cfl   - Execute the CFL check.  If the CFL check fails, then it will cause the program to fail. 
        The CFL check will check both the stability and accuracy conditions for p-waves. 
        
        NOTE: the CFL condition will return both minimum and maximum
        constraints on the grid given your velocity model, desired frequency content, and other
        parameters.  IT IS POSSIBLE TO HAVE NO STABLE, AND ACCURATE SOLUTIONS FOR A GIVEN 
        MODEL WITH GIVEN PARAMETERS. THE CFL CONDITION WILL WARN YOU IF THIS IS THE CASE.

        YOU MUST SPECIFY fmax Parameter as well!

         ----- STABILITY ------
        The stability condition is related to the maximum wave speed and minimum grid sampling
        as follows:

        dt < min(dx,dy,dz) / (sqrt(2)*vmax)

        Given a time sampling dt, it is possible to determine the minimum dx,dy,dz for stability.
        vmax is the MAXIMUM velocity for all waves in the model (usually P-wave velocity).

        The stability condition gives us a LOWER bound on the grid sampling for a given dt.

        ------ ACCURACY -------
        The accuracy condition is related to the number of gridpoints per wavelength.  Thus,

        safety*vmin / fmax > N * sqrt(dx^2+dy^2+dz^2) 

        where vmin is the minimum wave velocity in the model , fmax is some
        relative measure of the maximum frequency of your wavelet (usually 1.5*peak for Ricker), 
        N is the number of points desired per wavelength (5), and safety is a safety factor that 
        can be variable.

        The accuracy condition places an UPPER bound on the grid sampling.

        You can also override the safety factor using the safety parameter.
safety- Override the safety factor for the CFL condition.  This should be a floating point (0-1.0).
fmax  - An estimate of the highest frequency content in your wavelet (for Ricker use 1.5*peak)

fsrf  - Use a free surface at the top boundary (z=0).  

snap  - Save snapshots of the wavefield every n iterations of the modeling program. 

jsnap - Number of iterations between snapshots of the wavefield.  
	    i.e. jsnap=10, means save a snapshot every 10 iterations. 
	    If you had 1000 total iterations, then you would have 100 snapshots total.
	    The default, will output no snapshots.

jdata - Number of time imterations between exporting the data at the receivers.
	    i.e. jdata=1, means save a snapshot every iteration, which should be the default.
	    This can be used to change the sampling of the data to something different from 
        the wavelet/wavefield.

verb  - Print useful information

expl - Use an exploding reflector model (not sure what this does exactly as of 9/30/2010)

srctype - An integer which determines where the source wavelet is injected
        in the simulation.  Valid options are:  
            0 - Acceleration source
            1 - Displacement source
        The default option is 0: Acceleration source.

dabc  - Use a sponge layer to attenuate waves off the edge of the grid.  Use this in 
        combination with the nb parameter.

abcone- In addition to the sponge layer, using a severe ramp at the very edge of the expanded 
        sponge layer to severely attenuate zero-incidence waves at the boundaries. 
        It's not clear if this condition actually affects most computations.

nb    - Not listed, but is an important parameter.  Allows you to control the size of the sponge 
        layer for the absorbing boundary condition.  If you are getting reflections off the sides, 
        with dabc=y, then make this number larger (int).  This pads the grid by this amount on all sides.  
        For example:
                   |--------------------------|
                   |            ramp layer    |
                   |r |--------------------|  |
                   |a |        nb          |r |
                   |m |      |~~~~~~~~|    |a |
                   |p |      |  MODEL |    |m |
                   |  |  nb  |  SPACE | nb |p |
                   |  |      |~~~~~~~~|    |  |
                   |  |         nb         |  |
                   |  |--------------------|  |
                   |         ramp layer       |
                   |--------------------------| 
nqz, nqx, oqz, oqx, nqy, oqy, - Allows you to set the parameters for the axes.  Allows you to window the wavefield as it is output instead of having to do so afterwards..

====== BOUNDARY CONDITIONS =====

This program implements fixed, reflecting boundaries at the edge of the model space.


====== FILE DESCRIPTIONS ======

For filese that describe spatial values, the order of the spatial axes is ALWAYS:
    Axis 1- Z
    Axis 2 - X
    (Axis 3 - Y) (if necessary)

vp      - An N dimensional RSF file that contains the values for the velocity to be used for the model.  
           For 2D, this would be a 2D array.   

den      - An N dimensional RSF file that contains the valuese for the density to be used for the model.  
           For 2D, this would be a 2D array.  

sou, rec -The source and receiver RSF files respectively.  
          The 1st axis contains the locations for the points like so:
          [x,y,z]
          The second axis is a concatenated list of all points in the list.
          So, for an array of receivers, it would look like:
              [x1,y1,z1]
              [x2,y2,z2]
              [x3,y3,z3]
              [x4,y4,z4]

          In 2D, it will just be an array of [x,z] points.

wfl     - The name of the file to save the wavefield snapshots to.  This will be an N+2 
          dimensional file.  The file will be organized as follows:
              1-2(3) axes, spatial coordinates in z,x,y order
              3(4) axis, time, sequential snapshots
              ***The parentheses indicate what the axes will be for 3D models.

Fdat.rsf - An RSF file containing your data in the following format:
            axis 1 - source location
            axis 2 - Time

Fwav.rsf - An RSF file containing your wavelet information.  The wavelet needs 
           to have 3 samples on N1 one for each component Z-X-Y (or just Z-X for 2D).  The second 
           axis describes the component as a function of time.  The sampling interval, origin time, 
           and number of time samples will be used as the defaults for the modeling code.
	       i.e. your wavelet needs to have the same length and parameters that you want to model with!
	   Ex:
	   1st axis    index
	   Z component  0     0 0 0 0 0 0 0 0 1 2 3 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
	   X component  1     0 0 0 0 0 0 0 0 1 2 3 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
	   Y component  2     0 0 0 0 0 0 0 0 1 2 3 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
				    2nd axis

====== TROUBLESHOOTING ========


*/
/*
  Copyright (C) 2007 Colorado School of Mines
  
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
#include <stdio.h>
#include "fdutil.h"

/* check: dt<= 0.2 * min(dx,dz)/vmin */
#define NUM_INTERVALS 4 /* Number of intervals to use for CFL condition */
#define NOP 2 /* derivative operator half-size */

#define C0 -2.500000 /*    c0=-30./12.; */
#define CA +1.333333 /*    ca=+16./12.; */
#define CB -0.083333 /*    cb=- 1./12.; */

#define C1  0.66666666666666666666 /*  2/3  */	
#define C2 -0.08333333333333333333 /* -1/12 */

/* centered FD derivative stencils */
#define DX(a,ix,iz,s) (C2*(a[ix+2][iz  ] - a[ix-2][iz  ]) +  \
                       C1*(a[ix+1][iz  ] - a[ix-1][iz  ])  )*s
#define DZ(a,ix,iz,s) (C2*(a[ix  ][iz+2] - a[ix  ][iz-2]) +  \
                       C1*(a[ix  ][iz+1] - a[ix  ][iz-1])  )*s

#define DX3(a,ix,iy,iz,s) (C2*(a[iy  ][ix+2][iz  ] - a[iy  ][ix-2][iz  ]) +  \
                          C1*(a[iy  ][ix+1][iz  ] - a[iy  ][ix-1][iz  ])  )*s

#define DY3(a,ix,iy,iz,s) (C2*(a[iy+2][ix  ][iz  ] - a[iy-2][ix  ][iz  ]) +  \
                          C1*(a[iy+1][ix  ][iz  ] - a[iy-1][ix  ][iz  ])  )*s

#define DZ3(a,ix,iy,iz,s) (C2*(a[iy  ][ix  ][iz+2] - a[iy  ][ix  ][iz-2]) +  \
                          C1*(a[iy  ][ix  ][iz+1] - a[iy  ][ix  ][iz-1])  )*s

enum SourceType {ACCELERATION=0, DISPLACEMENT=1};

int main(int argc, char* argv[])
{
    bool verb,fsrf,snap,expl,dabc,abcone,is2D,cfl; 
    bool ignore_interpolation = false; /* ignore interpolation for receivers - makes code faster, but only works when receivers are on grid points */
    int  jsnap,ntsnap,jdata;
    float fmax, safety;
    enum SourceType srctype;
    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fvel=NULL; /* velocity  */
    sf_file Fden=NULL; /* density   */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */
/*set all y variables to be either zero or null to avoid compiler warnings
about being uninitialized */
    /* cube axes */
    sf_axis at,az,ax,ay=NULL; 
    sf_axis as,ar;

    int     nt,nz,nx,ny=0,ns,nr,nb;
    int     it,iz,ix,iy=0;
    float   dt,dz,dx,dy=0,idz,idx,idy=0;

    
    /* I/O arrays */
    float  *ww=NULL;           /* wavelet   */
    float  *dd=NULL;           /* data      */

    
    /* FD operator size */
    float co,cax,cbx,cay,cby,caz,cbz;

    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL,acy=NULL;
    int       nqz,nqx,nqy;
    float     oqz,oqx,oqy;
    float     dqz,dqx,dqy;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* OMP parameters */

    if( !sf_getbool("ignint",&ignore_interpolation)) ignore_interpolation = false;
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getbool("expl",&expl)) expl=false; /* "exploding reflector" */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* absorbing BC */
    if(! sf_getbool("cfl",&cfl)) cfl=false; /* Use CFL check */ 
    if(! sf_getbool("abcone",&abcone)) abcone=false; /* Use Zero-incident boundary condition*/ 
   
    int ttype = 0;
    if(! sf_getint("srctype",&ttype)) ttype = 0; /* source type, see comments */
    if(ttype < 0 || ttype > 1) sf_error("Invalid source type specified");
           srctype = ttype;
    if (cfl) {
        if(! sf_getfloat("fmax",&fmax)) { /* max frequency for cfl check */
            sf_error("CFL: Must specify fmax for CFL check");
        }
        if(! sf_getfloat("safety",&safety) || safety < 0.0) safety= 0.8; /*safety factor for cfl check*/
    }
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("in" ); /* wavelet   */
    Fvel = sf_input ("vel"); /* velocity  */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */
    Fden = sf_input ("den"); /* density   */

	/* Determine dimensionality, if 2D then axis 3 has n size of 1 */
	sf_axis test = sf_iaxa(Fvel,3);
	if(sf_n(test) == 1) is2D = true;
	else is2D = false;

    /*------------------------------------------------------------*/
    /* axes */
    at = sf_iaxa(Fwav,2); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    az = sf_iaxa(Fvel,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fvel,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */

    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */

    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);

    ns = sf_n(as);
    nr = sf_n(ar);

    if(!is2D){ /*If 3D*/
		ay=sf_iaxa(Fvel,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay); /*space*/
		ny=sf_n(ay); dy=sf_d(ay);
	}
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* other execution parameters */
    if(! sf_getint("jdata",&jdata)) jdata=1;
    if(snap) {  /* save wavefield every *jsnap* time steps */
    	if(! sf_getint("jsnap",&jsnap)) jsnap=nt;        
    }
    /*------------------------------------------------------------*/
if(is2D){
/* Begin 2d code */
    /* FDM structure */
    fdm2d    fdm=NULL;
    abcone2d abc=NULL;
    sponge   spo=NULL;
    pt2d   *ss=NULL;           /* sources   */
    pt2d   *rr=NULL;           /* receivers */
   
    float **tt=NULL;
    float **ro=NULL;           /* density */
    float **roz=NULL;          /* normalized 1st derivative of density on axis 1 */
    float **rox=NULL;          /* normalized 1st derivative of density on axis 2 */
    float **vp=NULL;           /* velocity */
    float **vt=NULL;           /* temporary vp*vp * dt*dt */

    float **um,**uo,**up,**ua,**ut; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */

    /* linear interpolation weights/indices */
    lint2d cs,cr;
    float     **uc=NULL;

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil_init(verb,fsrf,az,ax,nb,1);

    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* setup output data header */
    sf_oaxa(Fdat,ar,1);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,2);

    /* setup output wavefield header */
    if(snap) {
	if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az);
	if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax);

	if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
	if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
	dqz=sf_d(az);
	dqx=sf_d(ax);

	acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
	/* check if the imaging window fits in the wavefield domain */

	uc=sf_floatalloc2(sf_n(acz),sf_n(acx));

	ntsnap=0;
	for(it=0; it<nt; it++) {
	    if(it%jsnap==0) ntsnap++;
	}
	sf_setn(at,  ntsnap);
	sf_setd(at,dt*jsnap);
	if(verb) sf_raxa(at);

	sf_oaxa(Fwfl,acz,1);
	sf_oaxa(Fwfl,acx,2);
	sf_oaxa(Fwfl,at, 3);
    }

    if(expl) {
    	ww = sf_floatalloc( 1);
    } else {
    	ww = sf_floatalloc(ns);
    }
    dd = sf_floatalloc(nr);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt2d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(nr,sizeof(*rr)); 

    pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
    pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

    cs = lint2d_make(ns,ss,fdm);
    cr = lint2d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/
    /* setup FD coefficients */
    idz = 1/dz;
    idx = 1/dx;

    co = C0 * (idx*idx+idz*idz);
    cax= CA *  idx*idx;
    cbx= CB *  idx*idx;
    caz= CA *  idz*idz;
    cbz= CB *  idz*idz;

    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc2(nz,nx); 

    ro  =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    roz =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    rox =sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    vp  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
    vt  =sf_floatalloc2(fdm->nzpad,fdm->nxpad); 

    /* input density */
    sf_floatread(tt[0],nz*nx,Fden);     expand(tt,ro ,fdm);
    /* normalized density derivatives */
    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
	    roz[ix][iz] = DZ(ro,ix,iz,idz) / ro[ix][iz];
	    rox[ix][iz] = DX(ro,ix,iz,idx) / ro[ix][iz];
	}
    }   
    free(*ro); free(ro);

    /* input velocity */
    sf_floatread(tt[0],nz*nx,Fvel );    expand(tt,vp,fdm);
    float vpmax = 0.0; float vpmin = 10000000000000000;
    /* precompute vp^2 * dt^2 */
    for    (ix=0; ix<fdm->nxpad; ix++) {
        for(iz=0; iz<fdm->nzpad; iz++) {
            vt[ix][iz] = vp[ix][iz] * vp[ix][iz] * dt*dt;
            if (vp[ix][iz] < vpmin) vpmin = vp[ix][iz];
            else if (vp[ix][iz] > vpmax) vpmax = vp[ix][iz];
        }
    }
    if (cfl) cfl_acoustic(vpmin,vpmax,dx,-1.0f,dz,dt,fmax,safety,NUM_INTERVALS);
    if(fsrf) { /* free surface */
        for    (ix=0; ix<fdm->nxpad; ix++) {
            for(iz=0; iz<fdm->nb; iz++) {
                vt[ix][iz]=0;
            }
        }
    }

    free(*tt); free(tt);    
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    um=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uo=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    up=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    ua=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    for    (ix=0; ix<fdm->nxpad; ix++) {
	for(iz=0; iz<fdm->nzpad; iz++) {
	    um[ix][iz]=0;
	    uo[ix][iz]=0;
	    up[ix][iz]=0;
	    ua[ix][iz]=0;
	}
    }

    /*------------------------------------------------------------*/
	if (abcone) abc = abcone2d_make(NOP,dt,vp,fsrf,fdm);
    if(dabc) {
	/* one-way abc setup */
	/* sponge abc setup */
	spo = sponge_make(fdm->nb);
    }

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"%d/%d \r",it,nt);

#pragma omp parallel for				\
    schedule(dynamic) \
    private(ix,iz)					\
    shared(fdm,ua,uo,co,caz,cbz)
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		
		/* 4th order Laplacian operator */
		ua[ix][iz] = 
		    co * uo[ix  ][iz  ] + 
		    caz*(uo[ix  ][iz-1] + uo[ix  ][iz+1]) +
		    cbz*(uo[ix  ][iz-2] + uo[ix  ][iz+2]) ; 
		/* density term */
        /*ua[ix][iz] -= (
		    DZ(uo,ix,iz,idz) * roz[ix][iz] +
		    DX(uo,ix,iz,idx) * rox[ix][iz] );
        */
	    }
	}   

#pragma omp parallel for				\
    schedule(dynamic)			\
    private(ix,iz)					\
    shared(fdm,ua,uo,co,cax,cbx)
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
        ua[ix][iz] =  ua[ix][iz] + 
		    cax*(uo[ix-1][iz  ] + uo[ix+1][iz  ]) +
		    cbx*(uo[ix-2][iz  ] + uo[ix+2][iz  ]);
            }
    }

	/* inject acceleration source */
    if (srctype == ACCELERATION){
            if(expl) {
                sf_floatread(ww, 1,Fwav);
                lint2d_inject1(ua,ww[0],cs);
            } else {
                  sf_floatread(ww,ns,Fwav);
                  lint2d_inject(ua,ww,cs);
            }
    }

	/* step forward in time */
#pragma omp parallel for	    \
    schedule(dynamic) \
    private(ix,iz)		    \
    shared(fdm,ua,uo,um,up,vt)
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		up[ix][iz] = 2*uo[ix][iz] 
		    -          um[ix][iz] 
		    +          ua[ix][iz] * vt[ix][iz];
	    }
	}

    if(srctype == DISPLACEMENT){
            if(expl) {
                sf_floatread(ww, 1,Fwav);
                lint2d_inject1(up,ww[0],cs);
            } else {
                  sf_floatread(ww,ns,Fwav);
                  lint2d_inject(up,ww,cs);
            }
    }

    
	/* circulate wavefield arrays */
	ut=um;
	um=uo;
	uo=up;
	up=ut;
    
    if (abcone) abcone2d_apply(uo,um,NOP,abc,fdm);
	if(dabc) {
	    /* one-way abc apply */
	    sponge2d_apply(um,spo,fdm);
	    sponge2d_apply(uo,spo,fdm);
	    sponge2d_apply(up,spo,fdm);
	}

	/* extract data */
    if(ignore_interpolation){
        cut2d_extract(uo,dd,cr);
    } else {
	    lint2d_extract(uo,dd,cr);
    }

	if(snap && it%jsnap==0) {
	    cut2d(uo,uc,fdm,acz,acx);
	    sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
	}
	if(        it%jdata==0) 
	    sf_floatwrite(dd,nr,Fdat);
    }
    if(verb) fprintf(stderr,"\n");    

    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(*um); free(um);
    free(*up); free(up);
    free(*uo); free(uo);
    free(*ua); free(ua);
    if(snap) { free(*uc); free(uc); }

    free(*rox); free(rox);
    free(*roz); free(roz);
    free(*vp);  free(vp);
    free(*vt);  free(vt);

    free(ww);
    free(ss);
    free(rr);
    free(dd);


    exit (0);
} else {
    /* FDM structure */
    fdm3d    fdm=NULL;
    abcone3d abc=NULL;
    sponge   spo=NULL;

    /* I/O arrays */
    pt3d   *ss=NULL;           /* sources   */
    pt3d   *rr=NULL;           /* receivers */
	/* Non-universal arrays */
    float***tt=NULL;
    float***ro=NULL;           /* density */
    float***roz=NULL;          /* normalized 1st derivative of density on axis 1 */
    float***rox=NULL;          /* normalized 1st derivative of density on axis 2 */
    float***roy=NULL;          /* normalized 1st derivative of density on axis 3 */
    float***vp=NULL;           /* velocity */
    float***vt=NULL;           /* temporary vp*vp * dt*dt */

    float***um,***uo,***up,***ua,***ut; /* wavefield: um = U @ t-1; uo = U @ t; up = U @ t+1 */

    /* linear interpolation weights/indices */
    lint3d cs,cr;

	/* Wavefield cut params that are not universal */
    float     ***uc=NULL;


    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if(verb) sf_raxa(ay);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* setup output data header */
    sf_oaxa(Fdat,ar,1);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,2);

    /* setup output wavefield header */
    if(snap) {
	if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az);
	if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax);
	if(!sf_getint  ("nqy",&nqy)) nqy=sf_n(ay);

	if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
	if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);
	if(!sf_getfloat("oqy",&oqy)) oqy=sf_o(ay);

	dqz=sf_d(az);
	dqx=sf_d(ax);
	dqy=sf_d(ay);

	acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
	acy = sf_maxa(nqy,oqy,dqy); sf_raxa(acy);
	/* check if the imaging window fits in the wavefield domain */

	uc=sf_floatalloc3(sf_n(acz),sf_n(acx),sf_n(acy));

	ntsnap=0;
	for(it=0; it<nt; it++) {
	    if(it%jsnap==0) ntsnap++;
	}
	sf_setn(at,  ntsnap);
	sf_setd(at,dt*jsnap);
	if(verb) sf_raxa(at);

	sf_oaxa(Fwfl,acz,1);
	sf_oaxa(Fwfl,acx,2);
	sf_oaxa(Fwfl,acy,3);
	sf_oaxa(Fwfl,at, 4);
    }

    if(expl) {
	ww = sf_floatalloc( 1);
    } else {
	ww = sf_floatalloc(ns);
    }
    dd = sf_floatalloc(nr);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt3d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt3d*) sf_alloc(nr,sizeof(*rr)); 

    pt3dread1(Fsou,ss,ns,3); /* read (x,y,z) coordinates */
    pt3dread1(Frec,rr,nr,3); /* read (x,y,z) coordinates */

    cs = lint3d_make(ns,ss,fdm);
    cr = lint3d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/
    /* setup FD coefficients */
    idz = 1/dz;
    idx = 1/dx;
    idy = 1/dy;

    co = C0 * (idx*idx+idy*idy+idz*idz);
    cax= CA *  idx*idx;
    cbx= CB *  idx*idx;
    cay= CA *  idy*idy;
    cby= CB *  idy*idy;
    caz= CA *  idz*idz;
    cbz= CB *  idz*idz;

    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc3(nz,nx,ny); 

    ro  =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    roz =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    rox =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    roy =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    vp  =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
    vt  =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 

    /* input density */
    sf_floatread(tt[0][0],nz*nx*ny,Fden);     expand3d(tt,ro ,fdm);

    /* normalized density derivatives */
    for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		roz[iy][ix][iz] = DZ3(ro,ix,iy,iz,idz) / ro[iy][ix][iz];
		rox[iy][ix][iz] = DX3(ro,ix,iy,iz,idx) / ro[iy][ix][iz];
		roy[iy][ix][iz] = DY3(ro,ix,iy,iz,idy) / ro[iy][ix][iz];
	    }
	}   
    }
    free(**ro);  free(*ro); free(ro);  

    /* input velocity */
    sf_floatread(tt[0][0],nz*nx*ny,Fvel );    expand3d(tt,vp,fdm);
    /* precompute vp^2 * dt^2 */
    float vpmin = 1000000000000; float vpmax = 0.0;
    for        (iy=0; iy<fdm->nypad; iy++) {
        for    (ix=0; ix<fdm->nxpad; ix++) {
            for(iz=0; iz<fdm->nzpad; iz++) {
                float vpt = vp[iy][ix][iz];
                vt[iy][ix][iz] = vp[iy][ix][iz] * vp[iy][ix][iz] * dt*dt;
                if (vpt > vpmax) vpmax = vpt;
                else if (vpt < vpmin) vpmin = vpt;
            }
        }
    }

    if (cfl) cfl_acoustic(vpmin,vpmax,dx,dy,dz,dt,fmax,safety,NUM_INTERVALS);

    if(fsrf) { /* free surface */
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nb; iz++) {
		    vt[iy][ix][iz]=0;
		}
	    }
	}
    }

    free(**tt);  free(*tt); free(tt);    
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    um=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uo=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    up=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    ua=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    for        (iy=0; iy<fdm->nypad; iy++) {
	for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
		um[iy][ix][iz]=0;
		uo[iy][ix][iz]=0;
		up[iy][ix][iz]=0;
		ua[iy][ix][iz]=0;
	    }
	}
    }

    /*------------------------------------------------------------*/
	if (abcone) abc = abcone3d_make(NOP,dt,vp,fsrf,fdm);
    if(dabc) {
	/* one-way abc setup */
	/* sponge abc setup */
	spo = sponge_make(fdm->nb);
    }

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"%d/%d \r",it,nt);

#pragma omp parallel for					\
    schedule(dynamic) \
    private(ix,iy,iz)						\
    shared(fdm,ua,uo,co,cax,cay,caz,cbx,cby,cbz,idx,idy,idz)
	for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		    
		    /* 4th order Laplacian operator */
		    ua[iy][ix][iz] = 
			co * uo[iy  ][ix  ][iz  ] + 
			caz*(uo[iy  ][ix  ][iz-1] + uo[iy  ][ix  ][iz+1]) +
			cbz*(uo[iy  ][ix  ][iz-2] + uo[iy  ][ix  ][iz+2]);
		    
		    /* density term */
		    /*ua[iy][ix][iz] -= (
			DZ3(uo,ix,iy,iz,idz) * roz[iy][ix][iz] +
			DX3(uo,ix,iy,iz,idx) * rox[iy][ix][iz] +
			DY3(uo,ix,iy,iz,idy) * roy[iy][ix][iz] );
            */
		 }
	    }   
	}

#pragma omp parallel for					\
    schedule(dynamic) \
    private(ix,iy,iz)						\
    shared(fdm,ua,uo,co,cax,cay,caz,cbx,cby,cbz,idx,idy,idz)
	for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		    ua[iy][ix][iz] = ua[iy][ix][iz] + 
			cax*(uo[iy  ][ix-1][iz  ] + uo[iy  ][ix+1][iz  ]) +
			cbx*(uo[iy  ][ix-2][iz  ] + uo[iy  ][ix+2][iz  ]) ;
            }
        }
    }

#pragma omp parallel for					\
    schedule(dynamic) \
    private(ix,iy,iz)						\
    shared(fdm,ua,uo,co,cax,cay,caz,cbx,cby,cbz,idx,idy,idz)
	for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
		    ua[iy][ix][iz] = ua[iy][ix][iz] + 
			cay*(uo[iy-1][ix  ][iz  ] + uo[iy+1][ix  ][iz  ]) +
			cby*(uo[iy-2][ix  ][iz  ] + uo[iy+2][ix  ][iz  ]);
            }
        }
    }

	/* inject acceleration source */
    if (srctype == ACCELERATION){
        if(expl) {
            sf_floatread(ww, 1,Fwav);
            lint3d_inject1(ua,ww[0],cs);
        } else {
            sf_floatread(ww,ns,Fwav);	
            lint3d_inject(ua,ww,cs);
        }
   }

	/* step forward in time */
#pragma omp parallel for	    \
    schedule(static) \
    private(ix,iy,iz)		    \
    shared(fdm,ua,uo,um,up,vt)
	for        (iy=0; iy<fdm->nypad; iy++) {
	    for    (ix=0; ix<fdm->nxpad; ix++) {
		for(iz=0; iz<fdm->nzpad; iz++) {
		    up[iy][ix][iz] = 2*uo[iy][ix][iz] 
			-              um[iy][ix][iz] 
			+              ua[iy][ix][iz] * vt[iy][ix][iz];
		}
	    }
	}

    if (srctype == DISPLACEMENT) {
        if(expl) {
            sf_floatread(ww, 1,Fwav);
            lint3d_inject1(up,ww[0],cs);
        } else {
            sf_floatread(ww,ns,Fwav);	
            lint3d_inject(up,ww,cs);
        }
    }
	/* circulate wavefield arrays */
	ut=um;
	um=uo;
	uo=up;
	up=ut;

    if(abcone) abcone3d_apply(uo,um,NOP,abc,fdm);
	if(dabc) {
	    /* one-way abc apply */
	    sponge3d_apply(um,spo,fdm);
	    sponge3d_apply(uo,spo,fdm);
	    sponge3d_apply(up,spo,fdm);
	}

	/* extract data */
    if (ignore_interpolation) {
	    cut3d_extract(uo,dd,cr);
    } else {
	    lint3d_extract(uo,dd,cr);
    }

	if(snap && it%jsnap==0) {
	    cut3d(uo,uc,fdm,acz,acx,acy);
	    sf_floatwrite(uc[0][0],sf_n(acz)*sf_n(acx)*sf_n(acy),Fwfl);
	}
	if(it%jdata==0) 
	    sf_floatwrite(dd,nr,Fdat);
    }
    if(verb) fprintf(stderr,"\n");    

    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(**um); free(*um); free(um);
    free(**up); free(*up); free(up);
    free(**uo); free(*uo); free(uo);
    free(**ua); free(*ua); free(ua);
    if (snap) { free(**uc); free(*uc); free(uc); }

    free(**rox); free(*rox); free(rox);
    free(**roy); free(*roy); free(roy);
    free(**roz); free(*roz); free(roz);

    free(**vp); free(*vp); free(vp);
    free(**vt); free(*vt); free(vt);

    free(ss);
    free(rr);
    free(dd);
    free(ww);

    exit (0);
    }
}


