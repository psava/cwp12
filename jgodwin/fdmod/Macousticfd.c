#include <rsf.h>
#include "fdutility.h"
#include "fdacoustic.h"

/*



*/

int main(int argc, char* argv[]) 
{
    /* 
        Fire up RSF 
    */
    sf_init(argc,argv);
    /*
        Load all common FD modeling parameters into common
        datastructure for future use.  Add additional 
        parameters into this function.
    */
    fdutility_par par = fdutility_init_par(false);
    /*
        Initialize all the arrays that are needed for computation.
        Look in here for: array initialization and the CFL check.
    */
    fdacoustic_array array = fdacoustic_init_array(par);
    /*
        Fire off the simulation using the given parameters.
        For 2D, we dereference the 3D pointers for faster
        memory access and ignore the Y derivatives.
    */
    fdacoustic_simulate(array,par);
    /*
        Deallocate all of the memory that we have used.
    */
    fdacoustic_shutdown(array,par);
}
