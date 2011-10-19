#include <rsf.h>
#include "fdacoustic.h"
#include "fdutility.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// Type defs go here
#ifndef _fd_acoustic_h


typedef struct array *fdacoustic_array;

struct array {
    float ***uacce; //Acceleration 
    float ***ulast; //Last time step
    float ***uthis; //This time step
    float ***unext; //Next time step
    float ***usnap; //Snapshot wavefield (dims nqzXnqxXnqy)
    float ***rox;   // Density derivatives
    float ***roy;   // Density derivatives
    float ***roz;   // Density derivatives
    float ***vp;    // P-wave velocity
    float ***vt;    // Manipulated P-wave velocity
    float **ww;     // Wavelet 
};
/*^*/

#endif

void fdacoustic_simulate(fdacoustic_array array, fdutility_par par)
{

}



fdacoustic_array fdacoustic_init_arrays(fdacoustic_par par){
    fdacoustic_array array = (fdacoustic_array) sf_alloc(1,sizeof(*fdacoustic_array));

    
    if (par->verb) sf_warning("Allocating memory...\n");
    if (par->verb) sf_warning("Require: %d bytes \n",
        par->nxpad*par->nypad*par->nzpad*10*4);
    // Allocate the padded array
    float ***tt = sf_floatalloc3(par->nz,par->nx,par->ny);
    fdutility_zero_array(tt,par->nz,par->nx,par->ny);

    // Allocate the wavefield snapshot array
    array->usnap = sf_floatalloc3(par->nqz,par->nqx,par->nqy);
    fdutility_zero_array(array->usnap,par->nqz,par->nqx,par->nqy);

    // Allocate the model sized arrays and read them from file

    // Acceleration
    array->uacce = sf_floatalloc3(par->nzpad,par->nxpad,par->nypad);
    fdutility_zero_array(array->uacce,par->nzpad,par->nxpad,par->nypad);
   
    // Next time-step
    array->unext = sf_floatalloc3(par->nzpad,par->nxpad,par->nypad);
    fdutility_zero_array(array->unext,par->nzpad,par->nxpad,par->nypad);

    // This time-step
    array->uthis = sf_floatalloc3(par->nzpad,par->nxpad,par->nypad);
    fdutility_zero_array(array->uthis,par->nzpad,par->nxpad,par->nypad);
   
    // Last time-step
    array->ulast = sf_floatalloc3(par->nzpad,par->nxpad,par->nypad);
    fdutility_zero_array(array->ulast,par->nzpad,par->nxpad,par->nypad);

    // Velocity
    array->vp = sf_floatalloc3(par->nzpad,par->nxpad,par->nypad);
    fdutility_zero_array(array->vp,par->nzpad,par->nxpad,par->nypad);

    // Density derivative - X
    array->rox = sf_floatalloc3(par->nzpad,par->nxpad,par->nypad);
    fdutility_zero_array(array->rox,par->nzpad,par->nxpad,par->nypad);

    // Density derivative - Z
    array->roz = sf_floatalloc3(par->nzpad,par->nxpad,par->nypad);
    fdutility_zero_array(array->roz,par->nzpad,par->nxpad,par->nypad);

    if ( ! par->is2d) {
        // Density derivative - Y
        array->roy = sf_floatalloc3(par->nzpad,par->nxpad,par->nypad);
        fdutility_zero_array(array->roy,par->nzpad,par->nxpad,par->nypad);
    }

    /* 
    Store the density values in the VP array temporarily to save
    memory. We discard the ro values after we compute the gradient
    of the density anyways.
    */
    sf_floatread(tt[0][0],par->nx*par->ny*par->nz,par->Fden);
    fdutility_expand3d(tt,array->vp,par);

    for(int iy = 0; iy<par->ny; ++iy){
        for(int ix = 0; ix<par->nx; ++ix){
            for(int iz = 0; iz<par->nz; ++iz){
                rox[iy][ix][iz] = DX(vp,ix,iy,iz,idz)/vp[iy][ix][iz];
                roz[iy][ix][iz] = DZ(vp,ix,iy,iz,idz)/vp[iy][ix][iz];
                if (! par->is2d) roy[iy][ix][iz] = DY(vp,ix,iy,iz,idz)/vp[iy][ix][iz];
            }
        }
    }

    // Zero VP to make sure we don't keep anything from Density
    fdutility_zero_array(array->vp,par->nzpad,par->nxpad,par->nypad);
    // Now read the VP array into memory
    sf_floatread(tt[0][0],par->nx*par->ny*par->nz,par->Fvel);
    fdutility_expand3d(tt,array->vp,par);

    // Pre-compute scaled velocity for time-stepping
    // also grab the min and max velocity values in the 
    // model for use in the CFL check.
    float vpmin = 1e10; float vpmax = 0;
    for(int iy = 0; iy<par->ny; ++iy){
        for(int ix = 0; ix<par->nx; ++ix){
            for(int iz = 0; iz<par->nz; ++iz){
                if vp[iy][ix][iz] < vpmin: vpmin = vp[iy][ix][iz];
                if vp[iy][ix][iz] > vpmax: vpmax = vp[iy][ix][iz];
                vp[iy][ix][iz] = \
                    vp[iy][ix][iz]*vp[iy][ix][iz]*par->dt*par->dt;
            }
        }
    }
    // Read ALL wavelets into memory 
    // *** This could be large in 3D ***
    // Watch out here.
    array->ww = sf_floatalloc2(par->ns,par->nt);
    sf_floatread(ww,par->ns*par->nt,par->Fwav);

    /* 
     Now that we have the wavelets, find the frequency with maximum power
     in the frequency domain for all wavelet, and then find
     the value where the power is 1/2 the maximum for CFL check.
     i.e. for a Ricker wavelet with 20 Hz center frequency,
     find 40 Hz maximum frequency (where power is 1/2 peak).
    */
    int nw = par->nt/2+1;
    float dw = 1.0/(par->dt*par->nt);
    kiss_fft_cpx* cp = (kiss_fft_cpx*)sf_complexalloc(nw);
    kiss_fft_cfg cfg = sf_kiss_fft_alloc(par->nt,0,NULL,NULL);
    float *tw = sf_floatalloc(par->nt);
    for(int is = 0; is < par->ns; ++is){
        for(int it = 0; it < par->nt; ++it){
            tw[it] = ww[it][is];
        }
        kiss_fftr(cfg,tw,cp);
    }

    

}

void fdacoustic_shutdown(fdacoustic_array array)
{


}

