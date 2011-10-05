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
};
/*^*/

#endif

fdacoustic_array fdacoustic_init_arrays(fdacoustic_par par){
    fdacoustic_array array = (fdacoustic_array) sf_alloc(1,sizeof(*fdacoustic_array));

    // Allocate the padded array
    float ***tt = sf_floatalloc3(par->nz,par->nx,par->ny);
    fdutility_zero_array(tt,par->nz,par->nx,par->ny);

    // Allocate the wavefield snapshot array
    array->usnap = sf_floatalloc3(par->nqz,par->nqx,par->nqy);
    fdutility_zero_array(array->usnap,par->nqz,par->nqx,par->nqy);

    // Allocate the model sized arrays and read them from file
    array->uacce = sf_floatalloc3(par->nzpad,par->nxpad,par->nypad);
}

