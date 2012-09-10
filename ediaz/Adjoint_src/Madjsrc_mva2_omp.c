/* Computation of the adjoint source in the
   time domain using the L2 norm:

   J = || P(\tau) r(x,\tau)||^2_2

   then the adjoint source(s) should be:
   g_s = \sum_{x,\tau}r(x,\tau) P(\tau)^2 ur(x,t-2\tau)
   g_r = \sum_{x,\tau}r(x,\tau) P(\tau)^2 us(x,t+2\tau)




   source=0  apply a negative shift on the state variable u(x,t-2\tau)
   source=1  apply a positive shift on the state variable u(x,t+2\tau)



   This code only works with time-lags, but eventually
   will be extended for other lags.


   The 2-D wavefield (stdin) have the geometry from awefd2d (nz,nx,nt)

   The 2-D coordinates are read in duplets (x1,z1),(x2,z2)


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
#ifdef _OPENMP
#include <omp.h>
#endif



int main (int argc, char* argv[])
{
    /* ----------------------------------------------------------------------------------*/
    // Variable declaration
    // ------------------------
#ifdef _OPENMP
    int ompnth=1;
#endif


    // logical variables 
    bool source, verb;

    // axes
    sf_axis awx,awz,awt;                /* wavefield axes */
    sf_axis atau;                  /* eic axes: time-lag and extended images */
    sf_axis ar;                         /* coordinate file */

    int nr,nx,nz,nt;
    // integer
    int i1, i2,i3,it, itau,itaux;                   /* integer var for loops */

    int ix,iz;                          /*integers for extracting windowed wavefield */

    float tmin, tmax,taumin,taumax,t,tau,tbmin,tbmax;
    float *uaux;
    // arrays
    float ***uo;                        /* wavefield array */
    float *tgathers;                    /* time lag gathers */
    float **adjsrc;                      /* adjoint source [nt]  */

    // coordinate arrays
    pt2d *rr=NULL;                      /* extended images coordinates */

    /* ----------------------*/
    /* I/O files             */
    sf_file Feic=NULL;
    sf_file Fxcr=NULL;
    sf_file Fadj=NULL;
    sf_file Fwfl=NULL;
    /* ----------------------*/
    // End of variables declaration
    /* ----------------------------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /*------------------------------------------------------------*/
    /* OMP parameters */
#ifdef _OPENMP
    ompnth=omp_init();
    sf_warning("number of threads assigned by OMP: %d",ompnth);
#endif
    /*------------------------------------------------------------*/


    /* ------- I/O declaration  ---------------*/
    Fwfl = sf_input("in");      /* Wavefield used for the adjoint computation, for adj   
                                   source side pass the receiver wavefield, for receiver 
                                   side pass the source wavefield                        */
        
    Feic = sf_input("eic");     /* Penalized extended image (apply P(\tau) twice)        */

    Fxcr = sf_input("coord");   /* coordinates of every extended image, in 2d            
                                   should be organized in (x1,z1),(x2,z2),...,(xn,zn)    */

    Fadj = sf_output("out");    /* adjoint source to be injected on each coordinate      */

    /* --------------------------------------- */
    //Get parameters from command line:
    
    if (! sf_getbool("source",&source)) source=true; /* Source side [default] or receiver side adjoint?  */
    if (! sf_getbool("verb",&verb))       verb=true; /* verbosity flag, mainly debug printout */
                                                     /* at this time */
    
    if(verb) sf_warning("====================================================");
    if(verb) sf_warning("====== Calculating adjoint source            =======");
    if(source){ 
    if(verb) sf_warning("====== source side adjoint source            =======");
    }else{
    if(verb) sf_warning("====== receiver side adjoint source          =======");
    }
    if(verb) sf_warning("====================================================");
    // -------------------------------------------------
    // Some file type checking 
    if (SF_FLOAT != sf_gettype(Fwfl)) sf_error("Need float input for wavefield");
    if (SF_FLOAT != sf_gettype(Feic)) sf_error("Need float input for eic");
    if (SF_FLOAT != sf_gettype(Fxcr)) sf_error("Need float input for coordinates");
    // -------------------------------------------------

    // --------------------------------------------------------------------------------------
    //                      Axes geometry loading and writing
    // --------------------------------------------------------------------------------------
    // From wavefield
    awz = sf_iaxa(Fwfl,1); awx = sf_iaxa(Fwfl,2);  awt = sf_iaxa(Fwfl,3);

    // From coordinates
    ar = sf_iaxa(Fxcr,2); 
    
    // From extended images
    atau =  sf_iaxa(Feic,1);

    taumin =  SF_MIN(sf_o(atau),(sf_n(atau)-1)*sf_d(atau)+sf_o(atau));
    taumax =  SF_MAX(sf_o(atau),(sf_n(atau)-1)*sf_d(atau)+sf_o(atau));


    // To adjoint source
    sf_oaxa(Fadj,awt,1);
    sf_oaxa(Fadj,ar,2);
    sf_putint(Fadj,"n3",1);


    nr=sf_n(ar); 

    /*-------------------------------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    /*-------------------------------------------------------------------------------------*/
    rr = pt2dalloc1(nr);
    pt2dread1(Fxcr,rr,nr,2);




    if(verb) sf_warning("====================================================");
    if(verb) sf_warning("reading %d extended image coordinates points",nr);
    if(verb) sf_warning("====================================================");
    if(verb) sf_warning("allocating and reading wavefield");
    if(verb) sf_warning("z axis:  n1=%-10d o1=%-10.3f d1=%-10.3f",sf_n(awz),sf_o(awz),sf_d(awz));
    if(verb) sf_warning("x axis:  n2=%-10d o2=%-10.3f d2=%-10.3f",sf_n(awx),sf_o(awx),sf_d(awx));
    if(verb) sf_warning("t axis:  n3=%-10d o3=%-10.3f d3=%-10.3f",sf_n(awt),sf_o(awt),sf_d(awt));
    if(verb) sf_warning("====================================================");
   
    // allocate wavefield space    
    uo = sf_floatalloc3(sf_n(awz),sf_n(awx),sf_n(awt));

    if(uo==NULL) sf_error("not enough memory to read wavefield");

    uaux = sf_floatalloc(1);
    nz=sf_n(awz); nx=sf_n(awx); nt=sf_n(awt);

    // read wavefield
    for (i3=0 ; i3<nt; i3++){ 
        for (i2=0 ; i2<nx; i2++){
            for (i1=0 ; i1<nz; i1++){
                sf_floatread(uaux,1,Fwfl);
                uo[i3][i2][i1]=uaux[0];
    }}}




    //-------------------------------------------------------------------------------------
    // reading extended images 
    //-------------------------------------------------------------------------------------
    if(verb) sf_warning("reading %d extended images points",nr);
    if(verb) sf_warning("=====================================================");
    if(verb) sf_warning("");
    
    tgathers = sf_floatalloc(sf_n(atau));
    adjsrc   = sf_floatalloc2(sf_n(awt),1);



    if(source) {
        for (i1=0; i1<sf_n(ar); i1++) {
            sf_floatread(tgathers,sf_n(atau),Feic);

            ix=  0.5+(rr[i1].x - sf_o(awx))/sf_d(awx);
            iz=  0.5+(rr[i1].z - sf_o(awz))/sf_d(awz);




            fprintf(stderr,"\b\b\b\b\b\b%d",i1+1);
#ifdef _OPENMP
#pragma omp parallel for	    \
    private(it,t,tmin,tmax,itau,tbmin,tbmax,tau,itaux)		    \
    shared(tgathers,uo,)
#endif    
 
            for (it=0; it<sf_n(awt); it++){
    
                adjsrc[0][it]=0.0;

                // get the limits of tau such that we avoid
                // segmentation fault on the wavefield uo:
                // 0 < t-2*tau < tmax

                //Lower tau bound
                t=it*sf_d(awt)+sf_o(awt);

                tmin = SF_MIN(t+2*taumax,t+2*taumin);
                tmin = SF_MAX(tmin,0);
                 
                //Upper tau bound
        
                tmax = SF_MAX(t+2*taumax,t+2*taumin);
                tmax = SF_MIN(tmax,sf_o(awt)+sf_d(awt)*(sf_n(awt)-1) );


                tbmin =SF_MIN( -(t - tmax)*0.5,-(t - tmin)*0.5); 
                tbmax =SF_MAX( -(t - tmax)*0.5,-(t - tmin)*0.5);


                for (itau = 1.5 + (tbmin - sf_o(atau))/sf_d(atau)   ; itau<  (tbmax - sf_o(atau))/sf_d(atau)  ; itau++) {

                    tau=itau*sf_d(atau)+sf_o(atau) ;
                    itaux = 0.5 + ((t + 2*tau)-sf_o(awt))/sf_d(awt);

                    adjsrc[0][it] += tgathers[itau]*uo[itaux][ix][iz];
                }
            }
            sf_floatwrite(adjsrc[0],nt,Fadj);
        }  
    }else{

        for (i1=0; i1<sf_n(ar); i1++) {
            sf_floatread(tgathers,sf_n(atau),Feic);

            ix=  0.5+(rr[i1].x - sf_o(awx))/sf_d(awx);
            iz=  0.5+(rr[i1].z - sf_o(awz))/sf_d(awz);

            fprintf(stderr,"\b\b\b\b\b\b%d",i1+1);


            //if (verb) if(verb) fprintf(stderr,"\b\b\b\b\b%03d/%03d",i1+1,sf_n(ar));
#ifdef _OPENMP
#pragma omp parallel for	    \
    private(it,t,tmin,tmax,itau,tbmin,tbmax,tau,itaux)		    \
    shared(tgathers,uo,)
#endif    
            for (it=0; it<sf_n(awt); it++){
    
                adjsrc[0][it]=0.0;

                // get the limits of tau such that we avoid
                // segmentation fault on the wavefield uo:
                // 0 < t-2*tau < tmax

                //Lower tau bound
                t=it*sf_d(awt)+sf_o(awt);

                tmin = SF_MIN(t-2*taumax,t-2*taumin);
                tmin = SF_MAX(tmin,0);
                 
                //Upper tau bound
        
                tmax = SF_MAX(t-2*taumax,t-2*taumin);
                tmax = SF_MIN(tmax,sf_o(awt)+sf_d(awt)*(sf_n(awt)-1) );


                tbmin =SF_MIN( (t - tmax)*0.5,(t - tmin)*0.5); 
                tbmax =SF_MAX( (t - tmax)*0.5,(t - tmin)*0.5);


                for (itau = 1.5 + (tbmin - sf_o(atau))/sf_d(atau)   ; itau<  (tbmax - sf_o(atau))/sf_d(atau)  ; itau++) {

                    tau=itau*sf_d(atau)+sf_o(atau) ;
                    itaux = 0.5 + ((t - 2*tau)-sf_o(awt))/sf_d(awt);
                    
                    adjsrc[0][it] += tgathers[itau]*uo[itaux][ix][iz];
                }

            }
            sf_floatwrite(adjsrc[0],nt,Fadj);
        }  
    }

   if (verb) if(verb) fprintf(stderr,"\n"); 


    free(**uo);free(*uo);free(uo);
    free(*adjsrc);free(adjsrc);
    free(tgathers);
 
    exit(0);
}
