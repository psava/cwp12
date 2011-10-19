#include <rsf.h>
#include <stdio.h>
#include <math.h>
#include "fdutility.h"


#ifndef _fd_utility_h

typedef struct par  *fdutility_par;
/*^*/

typedef struct fdm3 *fdutility_model;
/*^*/

typedef struct abcone2 *fdutility_abcone2d;
/*^*/

typedef struct abcone3 *fdutility_abcone3d;
/*^*/

typedef struct sponge *fdutility_sponge;
/*^*/

typedef struct lcoef2 *fdutility_lint2d;
/*^*/

typedef struct lcoef3 *fdutility_lint3d;
/*^*/

struct par {
    bool verb;
    bool debug;
    bool snap;
    bool free;
    bool dabc;
    bool cfl;
    bool expl;
    bool abcone;
    int jsnap;
    int jdata;
    int srctype;
    float fmax;
    int nz, nx, ny, nt, ns, nr;
    int nzpad, nxpad, nypad;
    float oz, ox, oy, ot;
    float dz, dx, dy, dt;
    int nqz, nqx, nqy;
    float oqz, oqx, oqy;
    float dqz, dqx, dqy;
    int nb;
    int ani;
    sf_file Fwav;
    sf_file Fsou;
    sf_file Frec;
    sf_file Fvel;
    sf_file Fden;
    sf_file Fdat;
    sf_file Fwfl;
    bool is2d;
};
/*^*/

struct fdm3{
    int nb;
    int   nz,nzpad;
    int   nx,nxpad;
    int   ny,nypad;
    float oz,ozpad;
    float ox,oxpad;
    float oy,oypad;
    float dz;
    float dx;
    float dy;
};
/*^*/

struct lcoef2{
    int n;
    float *w00;
    float *w01;
    float *w10;
    float *w11;
    int *jz;
    int *jx;
    int nbell;
    float **bell;
};
/*^*/

struct lcoef3{
    int n;
    float *w000;
    float *w001;
    float *w010;
    float *w011;
    float *w100;
    float *w101;
    float *w110;
    float *w111;
    int *jz;
    int *jx;
    int *jy;
    int nbell;
    float ***bell;
};
/*^*/

struct abcone2{
    bool free;
    float *bzl;
    float *bzh;
    float *bxl;
    float *bxh;
};
/*^*/

struct abcone3{
    bool free;
    float**bzl;
    float**bzh;
    float**bxl;
    float**bxh;
    float**byl;
    float**byh;
};
/*^*/

struct sponge{
    float *w;
};
/*^*/

#endif

fdutility_model fdutility_model_init(int nb){
    
}


float*** fdutility_read_model(char *filename, fdutility_model model){
    sf_file Fmodel = sf_input(filename);
    int n1 = model->nz;
    int n2 = model->nx;
    int n3 = model->ny;

    float ***array = sf_floatalloc3(n1,n2,n3);
    sf_floatread(Fmodel,array,n1*n2*n3);
    sf_fileclose(Fmodel);
    return array;
}

void fdutililty_expand(float** a, 
	    float** b, 
	    fdm2d fdm)
/*< expand domain >*/
{
    int iz,ix;

    for     (ix=0;ix<fdm->nx;ix++) {
	for (iz=0;iz<fdm->nz;iz++) {
	    b[fdm->nb+ix][fdm->nb+iz] = a[ix][iz];
	}
    }

    for     (ix=0; ix<fdm->nxpad; ix++) {
	for (iz=0; iz<fdm->nb;    iz++) {
	    b[ix][           iz  ] = b[ix][           fdm->nb  ];
	    b[ix][fdm->nzpad-iz-1] = b[ix][fdm->nzpad-fdm->nb-1];
	}
    }

    for     (ix=0; ix<fdm->nb;    ix++) {
	for (iz=0; iz<fdm->nzpad; iz++) {
	    b[           ix  ][iz] = b[           fdm->nb  ][iz];
	    b[fdm->nxpad-ix-1][iz] = b[fdm->nxpad-fdm->nb-1][iz];
	}
    }
}

void fdutility_expand3d(float ***a, 
	      float ***b, 
	      fdutility_par fdm)
/*< expand domain >*/
{

    
    int iz,ix,i3;

    for         (i3=0;i3<fdm->ny;i3++) {
	for     (ix=0;ix<fdm->nx;ix++) {
	    for (iz=0;iz<fdm->nz;iz++) {
		b[fdm->nb+i3][fdm->nb+ix][fdm->nb+iz] = a[i3][ix][iz];
	    }
	}
    }

    for         (i3=0; i3<fdm->nypad; i3++) {
	for     (ix=0; ix<fdm->nxpad; ix++) {
	    for (iz=0; iz<fdm->nb;    iz++) {
		b[i3][ix][           iz  ] = b[i3][ix][           fdm->nb  ];
		b[i3][ix][fdm->nzpad-iz-1] = b[i3][ix][fdm->nzpad-fdm->nb-1];
	    }
	}
    }


    for         (i3=0; i3<fdm->nypad; i3++) {
	for     (ix=0; ix<fdm->nb;    ix++) {
	    for (iz=0; iz<fdm->nzpad; iz++) {
		b[i3][           ix  ][iz] = b[i3][           fdm->nb  ][iz];
		b[i3][fdm->nxpad-ix-1][iz] = b[i3][fdm->nxpad-fdm->nb-1][iz];
	    }
	}
    }

    for         (i3=0; i3<fdm->nb;    i3++) {
	for     (ix=0; ix<fdm->nxpad; ix++) {
	    for (iz=0; iz<fdm->nzpad; iz++) {
		b[           i3  ][ix][iz] = b[           fdm->nb  ][ix][iz];
		b[fdm->nypad-i3-1][ix][iz] = b[fdm->nypad-fdm->nb-1][ix][iz];
	    }
	}
    }

}


void fdutility_cut2d(float**  a,
	   float**  b,
	   fdm2d  fdm,
	   sf_axis c1, 
	   sf_axis c2)
/*< cut a rectangular wavefield subset >*/
{
    int iz,ix;
    int fz,fx;

    fz = (floor)((sf_o(c1)-fdm->ozpad)/fdm->dz);
    fx = (floor)((sf_o(c2)-fdm->oxpad)/fdm->dx);

    for     (ix=0;ix<sf_n(c2);ix++) {
	for (iz=0;iz<sf_n(c1);iz++) {
	    b[ix][iz] = a[fx+ix][fz+iz];
	}
    }
}

void fdutility_cut3d(float*** a,
	   float*** b,
	   fdm3d  fdm,
	   sf_axis c1, 
	   sf_axis c2,
	   sf_axis c3)
/*< cut a rectangular wavefield subset >*/
{
    int iz,ix,iy;
    int fz,fx,fy;

    fz = (floor)((sf_o(c1)-fdm->ozpad)/fdm->dz);
    fx = (floor)((sf_o(c2)-fdm->oxpad)/fdm->dx);
    fy = (floor)((sf_o(c3)-fdm->oypad)/fdm->dy);

    for         (iy=0;iy<sf_n(c3);iy++) {
	for     (ix=0;ix<sf_n(c2);ix++) {
	    for (iz=0;iz<sf_n(c1);iz++) {
		b[iy][ix][iz] = a[fy+iy][fx+ix][fz+iz];
	    }
	}
    }
}

/*------------------------------------------------------------*/
lint2d fdutility_lint2d_make(int    na, 
		   pt2d*  aa, 
		   fdm2d fdm)
/*< init 2D linear interpolation >*/
{
    lint2d ca;
    int    ia;
    float f1,f2;
    
    ca = (lint2d) sf_alloc(1,sizeof(*ca));

    ca->n = na;

    ca->w00 = sf_floatalloc(na);
    ca->w01 = sf_floatalloc(na);
    ca->w10 = sf_floatalloc(na);
    ca->w11 = sf_floatalloc(na);

    ca->jz  = sf_intalloc(na);
    ca->jx  = sf_intalloc(na);

    for (ia=0;ia<na;ia++) {
	
	if(aa[ia].z >= fdm->ozpad && 
	   aa[ia].z <  fdm->ozpad + (fdm->nzpad-1)*fdm->dz &&
	   aa[ia].x >= fdm->oxpad && 
	   aa[ia].x <  fdm->oxpad + (fdm->nxpad-1)*fdm->dx   ) {
	    
	    ca->jz[ia] = (int)( (aa[ia].z-fdm->ozpad)/fdm->dz);
	    ca->jx[ia] = (int)( (aa[ia].x-fdm->oxpad)/fdm->dx);
	    
	    f1 = (aa[ia].z-fdm->ozpad)/fdm->dz - ca->jz[ia];
	    f2 = (aa[ia].x-fdm->oxpad)/fdm->dx - ca->jx[ia];
	} else {
        sf_warning("YOUR SOURCES ARE OUTSIDE OF THE GRID!!!\n");
	    ca->jz[ia] = 0; 
	    ca->jx[ia] = 0;
	    
	    f1 = 1; 
	    f2 = 0;
	}

	ca->w00[ia] = (1-f1)*(1-f2);
	ca->w01[ia] = (  f1)*(1-f2);
	ca->w10[ia] = (1-f1)*(  f2);
	ca->w11[ia] = (  f1)*(  f2);
    }

    return ca;
}

/*------------------------------------------------------------*/
lint3d fdutility_lint3d_make(int    na, 
		   pt3d*  aa, 
		   fdm3d fdm)
/*< init 3D linear interpolation >*/
{
    lint3d ca;
    int    ia;
    float f1,f2,f3;
    
    ca = (lint3d) sf_alloc(1,sizeof(*ca));

    ca->n = na;

    ca->w000 = sf_floatalloc(na);
    ca->w001 = sf_floatalloc(na);
    ca->w010 = sf_floatalloc(na);
    ca->w011 = sf_floatalloc(na);
    ca->w100 = sf_floatalloc(na);
    ca->w101 = sf_floatalloc(na);
    ca->w110 = sf_floatalloc(na);
    ca->w111 = sf_floatalloc(na);

    ca->jz  = sf_intalloc(na);
    ca->jx  = sf_intalloc(na);
    ca->jy  = sf_intalloc(na);

    for (ia=0;ia<na;ia++) {
	
	if(aa[ia].z >= fdm->ozpad && 
	   aa[ia].z <  fdm->ozpad + (fdm->nzpad-1)*fdm->dz &&
	   aa[ia].x >= fdm->oxpad && 
	   aa[ia].x <  fdm->oxpad + (fdm->nxpad-1)*fdm->dx &&   
	   aa[ia].y >= fdm->oypad && 
	   aa[ia].y <  fdm->oypad + (fdm->nypad-1)*fdm->dy  ) {
	    
	    ca->jz[ia] = (int)( (aa[ia].z-fdm->ozpad)/fdm->dz);
	    ca->jx[ia] = (int)( (aa[ia].x-fdm->oxpad)/fdm->dx);
	    ca->jy[ia] = (int)( (aa[ia].y-fdm->oypad)/fdm->dy);
	    
	    f1 = (aa[ia].z-fdm->ozpad)/fdm->dz - ca->jz[ia];
	    f2 = (aa[ia].x-fdm->oxpad)/fdm->dx - ca->jx[ia];
	    f3 = (aa[ia].y-fdm->oypad)/fdm->dy - ca->jy[ia];

	} else {
	    ca->jz[ia] = 0; 
	    ca->jx[ia] = 0;
	    ca->jy[ia] = 0;
	    
	    f1 = 1; 
	    f2 = 0;
	    f3 = 0;
	}

	ca->w000[ia] = (1-f3)*(1-f1)*(1-f2);
	ca->w001[ia] = (1-f3)*(  f1)*(1-f2);
	ca->w010[ia] = (1-f3)*(1-f1)*(  f2);
	ca->w011[ia] = (1-f3)*(  f1)*(  f2);

	ca->w100[ia] = (  f3)*(1-f1)*(1-f2);
	ca->w101[ia] = (  f3)*(  f1)*(1-f2);
	ca->w110[ia] = (  f3)*(1-f1)*(  f2);
	ca->w111[ia] = (  f3)*(  f1)*(  f2);
    }

    return ca;
}


/*------------------------------------------------------------*/
void fdutility_lint2d_hold(float**uu,
		 float *ww,
		 lint2d ca)
/*< hold fixed value in field >*/
{
    int   ia;
    float wa;

#pragma omp parallel for schedule(dynamic,1) private(ia,wa) shared(ca,ww,uu)
    for (ia=0;ia<ca->n;ia++) {
	wa = ww[ia];
	
	uu[ ca->jx[ia]   ][ ca->jz[ia]   ] = wa;
	uu[ ca->jx[ia]   ][ ca->jz[ia]+1 ] = wa;
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]   ] = wa;
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] = wa;
    }
}

/*------------------------------------------------------------*/
void fdutility_lint2d_inject(float**uu,
		   float *ww,
		   lint2d ca)
/*< inject into wavefield >*/
{
    int   ia;
    float wa;

#pragma omp parallel for \
    schedule(dynamic,1) \
    private(ia,wa) \
    shared(ca,ww,uu)
    for (ia=0;ia<ca->n;ia++) {
	wa = ww[ia];
	
	uu[ ca->jx[ia]   ][ ca->jz[ia]   ] -= wa * ca->w00[ia];
	uu[ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= wa * ca->w01[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= wa * ca->w10[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= wa * ca->w11[ia];
    }
}

/*------------------------------------------------------------*/
void fdutility_lint3d_inject(float***uu,
		   float  *ww,
		   lint3d  ca)
/*< inject into wavefield >*/
{
    int   ia;
    float wa;

#pragma omp parallel for schedule(dynamic,1) private(ia,wa) shared(ca,ww,uu)
    for (ia=0;ia<ca->n;ia++) {
	wa = ww[ia];
	
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]   ] -= wa * ca->w000[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= wa * ca->w001[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= wa * ca->w010[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= wa * ca->w011[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]   ] -= wa * ca->w100[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= wa * ca->w101[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= wa * ca->w110[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= wa * ca->w111[ia];
    }
}

/*------------------------------------------------------------*/
void fdutility_lint2d_inject1(float**uu,
		    float  ww,
		    lint2d ca)
/*< inject into wavefield >*/
{
    int   ia;

#pragma omp parallel for schedule(dynamic,1) private(ia) shared(ca,ww,uu)
    for (ia=0;ia<ca->n;ia++) {

	uu[ ca->jx[ia]   ][ ca->jz[ia]   ] -= ww * ca->w00[ia];
	uu[ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= ww * ca->w01[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= ww * ca->w10[ia];
	uu[ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= ww * ca->w11[ia];
    }
}

/*------------------------------------------------------------*/
void fdutility_lint3d_inject1(float***uu,
		    float   ww,
		    lint3d  ca)
/*< inject into wavefield >*/
{
    int   ia;

#pragma omp parallel for schedule(dynamic,1) private(ia) shared(ca,ww,uu)
    for (ia=0;ia<ca->n;ia++) {

	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]   ] -= ww * ca->w000[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= ww * ca->w001[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= ww * ca->w010[ia];
	uu[ ca->jy[ia]   ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= ww * ca->w011[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]   ] -= ww * ca->w100[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]   ][ ca->jz[ia]+1 ] -= ww * ca->w101[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]   ] -= ww * ca->w110[ia];
	uu[ ca->jy[ia]+1 ][ ca->jx[ia]+1 ][ ca->jz[ia]+1 ] -= ww * ca->w111[ia];
    }
}

void fdutility_cut2d_extract(float **uu, float*dd, lint2d ca)
/*< extract from wavefield without interpolation>*/
{
    int ia;
//#pragma omp parallel for schedule(dynamic,1) private(ia) shared(ca,dd,uu)
    for (ia=0;ia<ca->n;ia++) {
	    dd[ia] = uu[ca->jx[ia]][ca->jz[ia]];
    }
}

/*------------------------------------------------------------*/
void fdutility_lint2d_extract(float**uu,
		    float* dd,
		    lint2d ca)
/*< extract from wavefield >*/
{
    int ia;

#pragma omp parallel for schedule(dynamic,1) private(ia) shared(ca,dd,uu)
    for (ia=0;ia<ca->n;ia++) {
	dd[ia] =
	    uu[ ca->jx[ia]  ][ ca->jz[ia]  ] * ca->w00[ia] +
	    uu[ ca->jx[ia]  ][ ca->jz[ia]+1] * ca->w01[ia] +
	    uu[ ca->jx[ia]+1][ ca->jz[ia]  ] * ca->w10[ia] +
	    uu[ ca->jx[ia]+1][ ca->jz[ia]+1] * ca->w11[ia];
    }
}  

void fdutility_lint3d_extract(float***uu,
		    float  *dd,
		    lint3d  ca)
/*< extract from wavefield >*/
{
    int ia;

#pragma omp parallel for schedule(dynamic,1) private(ia) shared(ca,dd,uu)
    for (ia=0;ia<ca->n;ia++) {
	dd[ia] =
	    uu[ ca->jy[ia]  ][ ca->jx[ia]  ][ ca->jz[ia]  ] * ca->w000[ia] +
	    uu[ ca->jy[ia]  ][ ca->jx[ia]  ][ ca->jz[ia]+1] * ca->w001[ia] +
	    uu[ ca->jy[ia]  ][ ca->jx[ia]+1][ ca->jz[ia]  ] * ca->w010[ia] +
	    uu[ ca->jy[ia]  ][ ca->jx[ia]+1][ ca->jz[ia]+1] * ca->w011[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]  ][ ca->jz[ia]  ] * ca->w100[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]  ][ ca->jz[ia]+1] * ca->w101[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]+1][ ca->jz[ia]  ] * ca->w110[ia] +
	    uu[ ca->jy[ia]+1][ ca->jx[ia]+1][ ca->jz[ia]+1] * ca->w111[ia];
    }
}  

void fdutility_cut3d_extract(float ***uu, float *dd, lint3d ca) {
/*< extract from wavefield without interpolation >*/
    int ia;
#pragma omp parallel for schedule(static) private(ia) shared(ca,dd,uu)
    for( ia=0; ia < ca->n; ++ia){
        dd[ia] = uu[ca->jy[ia]][ca->jx[ia]][ca->jz[ia]];
    }
}

/*------------------------------------------------------------*/
void fdutility_fdbell2d_init(int n)
/*< init bell taper >*/
{
    int   iz,ix;
    float s;

    nbell = n;
    s = 0.5*nbell;

    bell=sf_floatalloc2(2*nbell+1,2*nbell+1);

    for    (ix=-nbell;ix<=nbell;ix++) {
	for(iz=-nbell;iz<=nbell;iz++) {
	    bell[nbell+ix][nbell+iz] = exp(-(iz*iz+ix*ix)/s);
	}
    }    
}

/*------------------------------------------------------------*/
void fdutility_fdbell3d_init(int n)
/*< init bell taper >*/
{
    int   iz,ix,i3;
    float s;

    nbell = n;
    s = 0.5*nbell;

    bell3d=sf_floatalloc3(2*nbell+1,2*nbell+1,2*nbell+1);

    for        (i3=-nbell;i3<=nbell;i3++) {
	for    (ix=-nbell;ix<=nbell;ix++) {
	    for(iz=-nbell;iz<=nbell;iz++) {
		bell3d[nbell+i3][nbell+ix][nbell+iz] = exp(-(iz*iz+ix*ix+i3*i3)/s);
	    }
	}    
    }
}

/*------------------------------------------------------------*/
abcone2d fdutility_abcone2d_make(int     nop,
		       float    dt,
		       float**  vv,
		       bool   free, 
		       fdm2d   fdm)
/*< init 2D ABC >*/

/*This absorbing boundary condition follows the work done in
Absorbing Boundary Conditions , by Robert Clayton and Bjorn Engquist

Which can be found at:

http://sepwww.stanford.edu/public/docs/sep11/11_12_abs.html
*/
{
    abcone2d abc;
    int iz,ix;
    float d;

    abc = (abcone2d) sf_alloc(1,sizeof(*abc));

    abc->free = free;

    abc->bzl = sf_floatalloc(fdm->nxpad);
    abc->bzh = sf_floatalloc(fdm->nxpad);
    abc->bxl = sf_floatalloc(fdm->nzpad);
    abc->bxh = sf_floatalloc(fdm->nzpad);

    for (ix=0;ix<fdm->nxpad;ix++) {
	d = vv[ix][           nop  ] *dt/fdm->dz; abc->bzl[ix] = (1-d)/(1+d);
	d = vv[ix][fdm->nzpad-nop-1] *dt/fdm->dz; abc->bzh[ix] = (1-d)/(1+d);
    }
    for (iz=0;iz<fdm->nzpad;iz++) {
	d = vv[           nop  ][iz] *dt/fdm->dx; abc->bxl[iz] = (1-d)/(1+d);
	d = vv[fdm->nxpad-nop-1][iz] *dt/fdm->dx; abc->bxh[iz] = (1-d)/(1+d);
    }

    return abc;
}

/*------------------------------------------------------------*/
abcone3d fdutility_abcone3d_make(int     nop,
		       float    dt,
		       float ***vv,
		       bool   free, 
		       fdm3d   fdm)
/*< init 3D ABC >*/

/*This absorbing boundary condition follows the work done in
Absorbing Boundary Conditions , by Robert Clayton and Bjorn Engquist

Which can be found at:

http://sepwww.stanford.edu/public/docs/sep11/11_12_abs.html
*/

{
    abcone3d abc;
    int iz,ix,iy;
    float d;

    abc = (abcone3d) sf_alloc(1,sizeof(*abc));

    abc->free = free;

    /* z */
    abc->bzl = sf_floatalloc2(fdm->nxpad,fdm->nypad);
    abc->bzh = sf_floatalloc2(fdm->nxpad,fdm->nypad);

    /* x */
    abc->bxl = sf_floatalloc2(fdm->nzpad,fdm->nypad);
    abc->bxh = sf_floatalloc2(fdm->nzpad,fdm->nypad);

    /* y */
    abc->byl = sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    abc->byh = sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    for     (iy=0;iy<fdm->nypad;iy++) {
	for (ix=0;ix<fdm->nxpad;ix++) {
	    d = vv[iy][ix][           nop  ] *dt/fdm->dz; abc->bzl[iy][ix] = (1-d)/(1+d);
	    d = vv[iy][ix][fdm->nzpad-nop-1] *dt/fdm->dz; abc->bzh[iy][ix] = (1-d)/(1+d);
	}
    }

    for     (iy=0;iy<fdm->nypad;iy++) {
	for (iz=0;iz<fdm->nzpad;iz++) {
	    d = vv[iy][           nop  ][iz] *dt/fdm->dx; abc->bxl[iy][iz] = (1-d)/(1+d);
	    d = vv[iy][fdm->nxpad-nop-1][iz] *dt/fdm->dx; abc->bxh[iy][iz] = (1-d)/(1+d);
	}
    }
    
    for     (ix=0;ix<fdm->nxpad;ix++) {
	for (iz=0;iz<fdm->nzpad;iz++) {
	    d = vv[           nop  ][ix][iz] *dt/fdm->dy; abc->byl[ix][iz] = (1-d)/(1+d);
	    d = vv[fdm->nypad-nop-1][ix][iz] *dt/fdm->dy; abc->byh[ix][iz] = (1-d)/(1+d);
	}
    }

    return abc;
}


/*------------------------------------------------------------*/
void fdutility_abcone2d_apply(float**   uo,
		    float**   um,
		    int      nop,
		    abcone2d abc,
		    fdm2d    fdm)
/*< apply 2D ABC >*/
{
    int iz,ix,iop;

#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(iz,ix,iop)				\
    shared(fdm,nop,uo,um,abc)
    for(ix=0;ix<fdm->nxpad;ix++) {
	for(iop=0;iop<nop;iop++) {

	    /* top BC */
	    if(!abc->free) { /* not free surface, apply ABC */
		iz = nop-iop;
		uo      [ix][iz  ] 
		    = um[ix][iz+1] 
		    +(um[ix][iz  ]
		      - uo[ix][iz+1]) * abc->bzl[ix];
	    }

	    /* bottom BC */
	    iz = fdm->nzpad-nop+iop-1;
	    uo      [ix][iz  ] 
		= um[ix][iz-1]
		+(um[ix][iz  ]
		- uo[ix][iz-1]) * abc->bzh[ix];
	}
    }

#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(iz,ix,iop)				\
    shared(fdm,nop,uo,um,abc)
    for(iz=0;iz<fdm->nzpad;iz++) {
	for(iop=0;iop<nop;iop++) {

	    /* left BC */
	    ix = nop-iop;
	    uo      [ix  ][iz] 
		= um[ix+1][iz] 
		+(um[ix  ][iz]
		- uo[ix+1][iz]) * abc->bxl[iz];

	    /* right BC */
	    ix = fdm->nxpad-nop+iop-1;
	    uo      [ix  ][iz] 
		= um[ix-1][iz]
		+(um[ix  ][iz]
		- uo[ix-1][iz]) * abc->bxh[iz];
	}
    }
}

/*------------------------------------------------------------*/
void fdutility_abcone3d_apply(float  ***uo,
		    float  ***um,
		    int      nop,
		    abcone3d abc,
		    fdm3d    fdm)
/*< apply 3D ABC >*/
{
    int iz,ix,iy,iop;

#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(iz,ix,iy,iop)			\
    shared(fdm,nop,uo,um,abc)
    for    (iy=0;iy<fdm->nypad;iy++) {
	for(ix=0;ix<fdm->nxpad;ix++) {
	    for(iop=0;iop<nop;iop++) {
		
		/* z min */
		if(!abc->free) { /* not free surface, apply ABC */
		    iz = nop-iop;
		    uo      [iy][ix][iz  ] 
			= um[iy][ix][iz+1] 
			+(um[iy][ix][iz  ]
			- uo[iy][ix][iz+1]) * abc->bzl[iy][ix];
		}
		
		/* z max */
		iz = fdm->nzpad-nop+iop-1;
		uo      [iy][ix][iz  ] 
		    = um[iy][ix][iz-1]
		    +(um[iy][ix][iz  ]
		    - uo[iy][ix][iz-1]) * abc->bzh[iy][ix];
	    }
	}
    }
    
#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(iz,ix,iy,iop)			\
    shared(fdm,nop,uo,um,abc)
    for    (iy=0;iy<fdm->nypad;iy++) {
	for(iz=0;iz<fdm->nzpad;iz++) {
	    for(iop=0;iop<nop;iop++) {
		
		/* x min */
		ix = nop-iop;
		uo      [iy][ix  ][iz] 
		    = um[iy][ix+1][iz] 
		    +(um[iy][ix  ][iz]
		    - uo[iy][ix+1][iz]) * abc->bxl[iy][iz];
		
		/* x max */
		ix = fdm->nxpad-nop+iop-1;
		uo      [iy][ix  ][iz] 
		    = um[iy][ix-1][iz]
		    +(um[iy][ix  ][iz]
		    - uo[iy][ix-1][iz]) * abc->bxh[iy][iz];
	    }
	}
    }

#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(iz,ix,iy,iop)			\
    shared(fdm,nop,uo,um,abc)
    for    (ix=0;ix<fdm->nxpad;ix++) {
	for(iz=0;iz<fdm->nzpad;iz++) {
	    for(iop=0;iop<nop;iop++) {
		
		/* y min */
		iy = nop-iop;
		uo      [iy  ][ix][iz] 
		    = um[iy+1][ix][iz] 
		    +(um[iy  ][ix][iz]
		    - uo[iy+1][ix][iz]) * abc->byl[ix][iz];
		
		/* y max */
		iy = fdm->nypad-nop+iop-1;
		uo      [iy  ][ix][iz] 
		    = um[iy-1][ix][iz]
		    +(um[iy  ][ix][iz]
		    - uo[iy-1][ix][iz]) * abc->byh[ix][iz];
	    }
	}
    }

}


/*------------------------------------------------------------*/
sponge fdutility_sponge_make(int nb)
/*< init boundary sponge >*/

/* Sponge boundary conditions multiply incoming wavefields
by smaller coefficients to attenuate the wavefield over time and space.

The sponge coefficients need to deviate from 1 very gradually to ensure
that there are no induced reflections caused by large impedance 
contrasts */
{
    sponge spo;
    int   ib;
    float sb,fb;
    
    spo = (sponge) sf_alloc(1,sizeof(*spo));    
    spo->w = sf_floatalloc(nb);
    sb = 1.5*nb;               
    for(ib=0; ib<nb; ib++) {
	fb = ib/(sqrt(2.0)*sb);
	spo->w[ib] = exp(-fb*fb);
    }
    return spo;
}

/*------------------------------------------------------------*/
void fdutility_sponge2d_apply(float**   uu,
		    sponge   spo,
		    fdm2d    fdm)
/*< apply boundary sponge >*/
{
    int iz,ix,ib,ibz,ibx;
    float w;

//#pragma omp parallel for			\
    schedule(dynamic)				\
    private(ib,iz,ix,ibz,ibx,w)			\
    shared(fdm,uu)
    for(ix=0; ix<fdm->nxpad; ix++) {
        for(ib=0; ib<fdm->nb; ib++) {
        w = spo->w[fdm->nb-ib-1];
        uu[ix][ib ] *= w; /*    top sponge */
        }
    }

//#pragma omp parallel for			\
    schedule(dynamic)				\
    private(ib,iz,ix,ibz,ibx,w)			\
    shared(fdm,uu)
    for(ix=0; ix<fdm->nxpad; ix++) {
        for(ib=0; ib<fdm->nb; ib++) {
        w = spo->w[fdm->nb-ib-1];
        ibz = fdm->nzpad-ib-1;
        uu[ix][ibz] *= w; /* bottom sponge */
        }
    }

//#pragma omp parallel for			\
    schedule(dynamic)				\
    private(ib,iz,ix,ibz,ibx,w)			\
    shared(fdm,uu)
    for(ib=0; ib<fdm->nb; ib++) {
        w = spo->w[fdm->nb-ib-1];
        for(iz=0; iz<fdm->nzpad; iz++) {
            uu[ib ][iz] *= w; /*   left sponge */
        }
    }

//#pragma omp parallel for			\
    schedule(dynamic)				\
    private(ib,iz,ix,ibz,ibx,w)			\
    shared(fdm,uu)
    for(ib=0; ib<fdm->nb; ib++) {
        w = spo->w[fdm->nb-ib-1];
	    ibx = fdm->nxpad-ib-1;
        for(iz=0; iz<fdm->nzpad; iz++) {
            uu[ibx][iz] *= w; /*  right sponge */
        }
    }
}

void fdutility_sponge2d_apply_test(float**   uu,
		    sponge   spo,
		    fdm2d    fdm)
/*< apply boundary sponge >*/
{
    int iz,ix,ib,ibz,ibx;
    float w;

#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(ib,iz,ix,ibz,ibx,w)			\
    shared(fdm,uu)
    for(ib=0; ib<fdm->nb; ib++) {
	w = spo->w[ib];

	ibz = fdm->nzpad-ib-1;
	for(ix=0; ix<fdm->nxpad; ix++) {
	    uu[ix][ib ] *= w; /*    top sponge */
	    uu[ix][ibz] *= w; /* bottom sponge */
	}

	ibx = fdm->nxpad-ib-1;
	for(iz=0; iz<fdm->nzpad; iz++) {
	    uu[ib ][iz] *= w; /*   left sponge */
	    uu[ibx][iz] *= w; /*  right sponge */
	}

    }
}

/*------------------------------------------------------------*/
void fdutility_sponge3d_apply(float  ***uu,
		    sponge   spo,
		    fdm3d    fdm)
/*< apply boundary sponge >*/
{
    int iz,ix,iy,ib,ibz,ibx,iby;
    float w;

//#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(ib,iz,ix,iy,ibz,ibx,iby,w)		\
    shared(fdm,uu)
	for    (iy=0; iy<fdm->nypad; iy++) {
	    for(ix=0; ix<fdm->nxpad; ix++) {
            for(ib=0; ib<fdm->nb; ib++) {
                ibz = fdm->nzpad-ib-1;
                w = spo->w[fdm->nb-ib-1];
                uu[iy][ix][ib ] *= w; /* z min */
                uu[iy][ix][ibz] *= w; /* z max */
	        }
	    }
    }

//#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(ib,iz,ix,iy,ibz,ibx,iby,w)		\
    shared(fdm,uu)
	for    (iy=0; iy<fdm->nypad; iy++) {
        for(ib=0; ib<fdm->nb; ib++) {
	        ibx = fdm->nxpad-ib-1;
	        for(iz=0; iz<fdm->nzpad; iz++) {
		        uu[iy][ib ][iz] *= w; /* x min */
		        uu[iy][ibx][iz] *= w; /* x max */
            }
        }
    }
	
//#pragma omp parallel for			\
    schedule(dynamic,1)				\
    private(ib,iz,ix,iy,ibz,ibx,iby,w)		\
    shared(fdm,uu)
    for(ib=0; ib<fdm->nb; ib++) {
	    iby = fdm->nypad-ib-1;
        for    (ix=0; ix<fdm->nxpad; ix++) {
            for(iz=0; iz<fdm->nzpad; iz++) {
            uu[ib ][ix][iz] *= w; /* x min */
            uu[iby][ix][iz] *= w; /* x max */
            }
        }
    }
} 


/*------------------------------------------------------------*/
void fdutility_lint2d_bell(float**uu,
		 float *ww,
		 lint2d ca)
/*< apply bell taper >*/
{
    int   ia,iz,ix;
    float wa;

    for    (ix=-nbell;ix<=nbell;ix++) {
	for(iz=-nbell;iz<=nbell;iz++) {
	    
	    for (ia=0;ia<ca->n;ia++) {
		wa = ww[ia] * bell[nbell+ix][nbell+iz];

		uu[ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] -= wa * ca->w00[ia];
		uu[ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] -= wa * ca->w01[ia];
		uu[ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] -= wa * ca->w10[ia];
		uu[ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] -= wa * ca->w11[ia];
	    }

	}
    }
}

/*------------------------------------------------------------*/
void fdutility_lint3d_bell(float***uu,
		 float  *ww,
		 lint3d  ca)
/*< apply bell taper >*/
{
    int   ia,iz,ix,iy;
    float wa;

    for        (iy=-nbell;iy<=nbell;iy++) {
	for    (ix=-nbell;ix<=nbell;ix++) {
	    for(iz=-nbell;iz<=nbell;iz++) {
		
		for (ia=0;ia<ca->n;ia++) {
		    wa = ww[ia] * bell3d[nbell+iy][nbell+ix][nbell+iz];
		    
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] -= wa * ca->w000[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] -= wa * ca->w001[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] -= wa * ca->w010[ia];
		    uu[ iy+ca->jy[ia]   ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] -= wa * ca->w011[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]   ] -= wa * ca->w100[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]   ][ iz+ca->jz[ia]+1 ] -= wa * ca->w101[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]   ] -= wa * ca->w110[ia];
		    uu[ iy+ca->jy[ia]+1 ][ ix+ca->jx[ia]+1 ][ iz+ca->jz[ia]+1 ] -= wa * ca->w111[ia];
		}
		
	    }
	}
    }
}

fdutility_par fdutility_init_par(bool elastic)
{
    fdutility_par fdpar = (fdutility_par) sf_alloc(1,sizeof(*fdutility_par));
    if (! sf_getbool("verb",&fdpar->verb)) fdpar->verb = false; /* verbose output */
    if (! sf_getbool("debug",&fdpar->verb)) fdpar->debug= false; /* verbose output */
    if (! sf_getbool("snap",&fdpar->snap)) fdpar->snap = false; /* output wfld snapshots */
    if (! sf_getbool("free",&fdpar->free)) fdpar->free = false; /* free surface */
    if (! sf_getbool("expl",&fdpar->expl)) fdpar->expl = false; /* exploding reflector */
    if (! sf_getbool("cfl",&fdpar->cfl)) fdpar->cfl = false; /* use CFL and stability check?*/
    if (! sf_getbool("dabc",&fdpar->dabc)) fdpar->dabc = false; /* sponge boundary? */
    if (! sf_getbool("abcone",&fdpar->abcone)) fdpar->abcone = false; /* plane-wave destruction boundary?*/

    if (! sf_getint("jsnap",&fdpar->jsnap)) fdpar->jsnap = 1; /* output wfld every jsnap time steps */
    if (! sf_getint("jdata",&fdpar->jdata)) fdpar->jdata = 1; /* output data every jdata time steps */
    if (! sf_getint("srctype",&fdpar->srctype)) fdpar->srctype = 0; /* source type: 0 - acceleration, 1 - displacement*/
    if (! sf_getint("nb",&fdpar->nb)) fdpar->nb = 0; /* # of boundary cells on each side for sponge*/
    if (! sf_getint("ani",&fdpar->ani)) fdpar->ani = -1; /* Type of anisotropy (elastic only) */

    if(fdpar->ani == -1 && elastic){
        sf_error("Must specify type of anisotropy! ani=-1 right now\n");
    } else if (elastic) {
        sf_warning("Specified: %d anisotropy\n",ani);
    }

    if (fdpar->verb){
        sf_warning("Opening files...\n");
    }

    fdpar->Fwav = sf_input("in");
    fdpar->Fsou = sf_input("sou");
    fdpar->Frec = sf_input("rec");
    if (elastic) fdpar->Fvel = sf_input("ccc");
    else fdpar->Fvel = sf_input("vel");
    fdpar->Fden = sf_input("den");
    fdpar->Fdat = sf_output("out");
    if (fdpar->snap) fdpar->Fwav=sf_output("wfl");
    else fdpar->Fwav = NULL;

    sf_axis at = sf_iaxa(fdpar->Fwav,1); sf_setlabel(at,"t");
    sf_axis az = sf_iaxa(fdpar->Fvel,1); sf_setlabel(az,"z");
    sf_axis ax = sf_iaxa(fdpar->Fvel,2); sf_setlabel(ax,"x");
    sf_axis ay = sf_iaxa(fdpar->Fvel,3); sf_setlabel(ay,"y"); 
    sf_axis as = sf_iaxa(fdpar->Fsou,2); sf_setlabel(as,"sou");
    sf_axis ar = sf_iaxa(fdpar->Frec,2); sf_setlabel(ar,"rec");

    fdpar->nt = sf_n(at); fdpar->ot = sf_o(at); fdpar->dt = sf_d(at);
    fdpar->nz = sf_n(az); fdpar->oz = sf_o(az); fdpar->dz = sf_d(az);
    fdpar->nx = sf_n(ax); fdpar->ox = sf_o(ax); fdpar->dx = sf_d(ax);
    fdpar->ny = sf_n(ay); fdpar->oy = sf_o(ay); fdpar->dy = sf_d(ay);

    if (! sf_getint("nqz",&fdpar->nqz)) fdpar->nqz = sf_n(az); /* # of output samples in wfld snapshots */
    if (! sf_getint("nqx",&fdpar->nqx)) fdpar->nqx = sf_n(ax); /* # of output samples in wfld snapshots */
    if (! sf_getint("nqy",&fdpar->nqy)) fdpar->nqy = sf_n(ay); /* # of output samples in wfld snapshots */
    if (! sf_getfloat("dqz",&fdpar->dqz)) fdpar->dqz = sf_d(az); /* sampling interval for wfld snapshots */
    if (! sf_getfloat("dqx",&fdpar->dqx)) fdpar->dqx = sf_d(ax); /* sampling interval for wfld snapshots */
    if (! sf_getfloat("dqy",&fdpar->dqy)) fdpar->dqy = sf_d(ay); /* sampling interval for wfld snapshots */
    if (! sf_getfloat("oqz",&fdpar->oqz)) fdpar->oqz = sf_o(az); /* origin for wfld snapshots */
    if (! sf_getfloat("oqx",&fdpar->oqx)) fdpar->oqx = sf_o(ax); /* origin for wfld snapshots */
    if (! sf_getfloat("oqy",&fdpar->oqy)) fdpar->oqy = sf_o(ay); /* origin for wfld snapshots */

    if (fdpar->ny == 1) fdpar->is2d = true;
    else fdpar-> d3 = false;
    fdpar->nzpad = fdpar->nz+2*nb; //+2 is for the abcone boundary
    fdpar->nxpad = fdpar->nx+2*nb;
    if (! fdpar->is2d) fdpar->nypad = fdpar->ny+2*nb; 
    else fdpar->nypad = 1; //if 2D don't pad 
    if (fdpar->abcone) { //Pad with two extra cells so that we can run abcone at the boundary
        fdpar->nzpad += 2;
        fdpar->nxpad += 2;
        if (! fdpar->is2d) fdpar->nypad += 2;
    }



    sf_axis tdat = sf_maxa(fdpar->nt/fdpar->jdata,fdpar->ot,fdpar->dt*fdpar->jdata); sf_setlabel(tdat,"data t");
    sf_axis zwfl = sf_maxa(fdpar->nqz,fdpar->oqz,fdpar->dqz); sf_setlabel("z wfl");
    sf_axis xwfl = sf_maxa(fdpar->nqx,fdpar->oqx,fdpar->dqx); sf_setlabel("x wfl");
    sf_ayis ywfl = sf_maya(fdpar->nqy,fdpar->oqy,fdpar->dqy); sf_setlabel("y wfl");

    if(fdpar->verb){
        sf_warning("Input axes information: \n");
        sf_raxa(at);
        sf_raxa(az);
        sf_raxa(ax);
        sf_raxa(ay);
        sf_raxa(as);
        sf_raxa(ar);
        sf_warning("Output axes information: \n");
        sf_raxa(tdat);

        if (fdpar->snap){
            sf_raxa(zwfl);
            sf_raxa(xwfl);
            sf_raxa(ywfl);
        }
    }
    if(fdpar->snap){
        if(fdpar->verb) sf_warning("Checking snapshot coordinates...");
        float zmin = fdpar->oz; 
        float zmax = fdpar->oz+(fdpar->nz-1)*fdpar->dz;
        float xmin = fdpar->ox;
        float xmax = fdpar->ox+(fdpar->nx-1)*fdpar->dx;
        float ymin = fdpar->oy;
        float ymay = fdpar->oy+(fdpar->ny-1)*fdpar->dy;

        float qzmin = fdpar->oqz; 
        float qzmax = fdpar->oqz+(fdpar->nqz-1)*fdpar->dqz;
        float qxmin = fdpar->oqx;
        float qxmax = fdpar->oqx+(fdpar->nqx-1)*fdpar->dqx;
        float qymin = fdpar->oqy;
        float qymax = fdpar->oqy+(fdpar->nqy-1)*fdpar->dqy;

        if (qzmin >= zmin && qzmax <= zmax && \
           qxmin >= xmin && qxmax <= xmax && \
           qymin >= ymin && qymax <= ymax){
            if(fdpar->verb) sf_warning("Passed\n");
        } else {
            sf_error("Wavefield snapshot domain outside of wavefield computation domain\n");
        }
    }
    return fdpar;
}

void fdutility_zero_array(float ***array,int n1, int n2, int n3){
    for(int i3=0; i3 < n3; ++i3){
        for(int i2=0; i2 < n2; ++i2){
            for(int i1=0; i1 < n1; ++i1){
                array[i3][i2][i1] = 0;
            }
        }
    }
}
