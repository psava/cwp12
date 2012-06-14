/* Computation of the adjoint source in the
   time domain using the L2 norm:

   J = || P(\tau) r(x,\tau)||^2_2

   then the adjoint should be:
   g_s = \sum_{x,\tau}r(x,\tau) P(\tau) ur(x,t+\tau) 

   This code only works with time-lags, but eventually
   will be extended for other lags.

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
    /* ----------------------------------------------------------------------------------*/
    // Variable declaration
    // ------------------------

    // logical variables 
    bool source, verb;

    // axes
    sf_axis awx,awz,awt;                /* wavefield axes */
    sf_axis atau;                  /* eic axes: time-lag and extended images */
    sf_axis ar;                         /* coordinate file */

    int nr;
    // integer
    int i1, it, itau;                   /* integer var for loops */

    int ix,iz;                          /*integers for extracting windowed wavefield */

    float tmin, tmax;
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

    /* ----------------------*/
    // init RSF argument parsing
    sf_init (argc,argv);
    // end of init RSF argument parsing

    /* ------- I/O declaration  ---------------*/
    Fwfl = sf_input("in");      /* Wavefield used for the adjoint computation, for adj   */
                                /* source side pass the receiver wavefield, for receiver */
                                /* side pass the source wavefield                        */
        
    Feic = sf_input("eic");     /* Penalized extended image (apply P(\tau) twice)        */

    Fxcr = sf_input("coord");   /* coordinates of every extended image, in 2d            */
                                /* should be organized in (x1,z1),(x2,z2),...,(xn,zn)    */

    Fadj = sf_output("out");    /* adjoint source to be injected on each coordinate      */

    /* --------------------------------------- */
    //Get parameters from command line:
    
    if (! sf_getbool("source",&source)) source=true; /* Source side [default] or receiver side adjoint?  */
    if (! sf_getbool("verb",&verb))       verb=true; /* verbosity flag, mainly debug printout */
                                                     /* at this time */
    
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

    // To adjoint source
    sf_oaxa(Fadj,awt,1);
    sf_oaxa(Fadj,ar,2);
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

    // read wavefield
    sf_floatread(uo[0][0],sf_n(awz)*sf_n(awx)*sf_n(awt),Fwfl);

    /*-------------------------------------------------------------------------------------*/
    /* reading extended images */
    /*-------------------------------------------------------------------------------------*/
    if(verb) sf_warning("reading %d extended images points",nr);
    if(verb) sf_warning("=====================================================");
    if(verb) sf_warning("");
    atau =  sf_iaxa(Feic,1);
    
    tgathers = sf_floatalloc(sf_n(atau));
    adjsrc   = sf_floatalloc2(sf_n(awt),1);



    if(source) {
        for (i1=0; i1<sf_n(ar); i1++) {
            sf_floatread(tgathers,sf_n(atau),Feic);
    
            ix=  0.5+(rr[i1].x - sf_o(awx))/sf_d(awx);
            iz=  0.5+(rr[i1].z - sf_o(awz))/sf_d(awz);
    
            sf_warning("cc=%d %d",i1,sf_n(atau));
    
            for (it=0; it<sf_n(awt); it++){
    
            sf_warning("it=%d %d",it,sf_n(awt));
                adjsrc[0][it]=0.0;

                tmin=SF_MAX(sf_o(awt)+it*sf_d(awt)+sf_o(atau),0);
                tmax=SF_MIN(sf_o(awt)+it*sf_d(awt)+(sf_o(atau)+(sf_n(atau)-1)*sf_d(atau)),(sf_n(awt)-1)*sf_d(awt)+sf_o(awt));
                

                for (itau=0.5 + (tmin - sf_o(atau))/sf_d(atau)  ; itau<0.5 + (tmax - sf_o(atau))/sf_d(atau) ; itau++) {
                    sf_warning("                itau=%d %f",itau,sf_o(atau)+(itau)*sf_d(atau));
                    adjsrc[0][it] += tgathers[itau]*uo[it+itau][ix][iz];
                }
            }
        }   
    }else{
        for (i1=0; i1<sf_n(ar); i1++) {
            sf_floatread(tgathers,sf_n(atau),Feic);
    
            ix=  0.5+(rr[i1].x - sf_o(awx))/sf_d(awx);
            iz=  0.5+(rr[i1].z - sf_o(awz))/sf_d(awz);
    
            sf_warning("cc=%d %d",i1,sf_n(atau));

            for (it=0; it<sf_n(awt); it++){
    
                adjsrc[0][it]=0.0;

                tmin=SF_MAX(sf_o(awt)+it*sf_d(awt)-(sf_o(atau)+(sf_n(atau)-1)*sf_d(atau)),0);
                tmax=SF_MIN(sf_o(awt)+it*sf_d(awt)-(sf_o(atau)),(sf_n(awt)-1)*sf_d(awt)+sf_o(awt));
                

                for (itau=0.5 + (tmin - sf_o(atau))/sf_d(atau)  ; itau<0.5 + (tmax - sf_o(atau))/sf_d(atau) ; itau++) {
                    adjsrc[0][it] += tgathers[itau]*uo[it-itau][ix][iz];
                }
            }





        }   


    }
    exit(0);
}
