#include <rsf.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "omp.h"


int main(){
    int nthreads = 4;
    omp_set_num_threads(nthreads);
    #pragma omp parallel 
        fprintf(stderr,"nthreads %d \n", omp_get_num_threads());
 
    int n3 = 128;
    int n2 = 128;
    int n1 = 128;
//    float ***array = sf_floatalloc3(n1,n2,n3);
    
    float *array = fftwf_alloc_real(n3*n2*n1);
    fftwf_complex* cout = fftwf_alloc_complex(n3*n2*n1);

    int err = fftwf_init_threads();
    if (err == 0) {
        fprintf(stderr,"something went wrong with fftw\n");
    }

    fprintf(stderr,"Got here\n");

    double start,end;
    start = omp_get_wtime()*omp_get_wtick();
    fftwf_plan_with_nthreads(nthreads);
    fftwf_plan plan =  fftwf_plan_dft_r2c_3d(
                                    n1,n2,n3,
                                    array,cout,
                                    FFTW_MEASURE);
    end = omp_get_wtime()*omp_get_wtick();
    fprintf(stderr,"elapsed time: %f %f %f\n",end,start,end-start);

    for(int i = 0; i < n3*n2*n1; ++i)
        array[i] = rand()/RAND_MAX;
 
    //float start = clock()/CLOCKS_PER_SEC;
    start = omp_get_wtime();

    for(int i=0; i < 1001; ++i)
        fftwf_execute(plan);
   
    //float end = clock()/CLOCKS_PER_SEC;
    end = omp_get_wtime();
    fprintf(stderr,"elapsed time: %f time/calc %f\n",
        end-start,(end-start)/100.0);

    fftwf_cleanup_threads();
    fftwf_cleanup();
    fftwf_destroy_plan(plan);

    fftwf_free(cout);
    fftwf_free(array);
    //free(**array); free(*array); free(array);
    return 0;

}
