/*

  Extract the value of an image (input) at sparse 
  locations provided by the coordinate file (coord= ). 

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


    // axes
    sf_axis awx,awz,axcr;
    sf_axis at; 

    int nr,nx,nz;
    float oz,ox, dz,dx;

    // integer
    int i1;                   /* integer var for loops */

    int ix,iz;                       /*integers for extracting windowed wavefield */

    // arrays
    float **image;                   /* input image */
    float *sparse_image; 

    // coordinate arrays
    pt2d *rr=NULL;                      /* extended images coordinates */

    /* ----------------------*/
    /* I/O files             */
    sf_file Fimg=NULL;
    sf_file Fxcr=NULL;
    sf_file Fout=NULL;
    /* ----------------------*/
    // End of variables declaration
    /* ----------------------------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);



    /* ------- I/O declaration  ---------------*/
    Fimg = sf_input("in");      /* input image           */

    Fxcr = sf_input("coord");   /* coordinate file with sparse locations to extract  */

    Fout = sf_output("out");    /* output vector with sparse locations      */
    
    // -------------------------------------------------
    // Some file type checking 
    if (SF_FLOAT != sf_gettype(Fimg)) sf_error("Need float input for image");
    if (SF_FLOAT != sf_gettype(Fxcr)) sf_error("Need float input for coordinates");
    // -------------------------------------------------

    // --------------------------------------------------------------------------------------
    //                      Axes geometry loading and writing
    // --------------------------------------------------------------------------------------
    // From wavefield
    awz = sf_iaxa(Fimg,1); awx = sf_iaxa(Fimg,2);  

    nz = sf_n(awz); oz = sf_o(awz) ; dz = sf_d(awz);
    nx = sf_n(awx); ox = sf_o(awx) ; dx = sf_d(awx);


    // From coordinates
    axcr = sf_iaxa(Fxcr,2); 
    

    //=========================== 
    // write output header:
    at = sf_iaxa(Fxcr,1) ;
    
    sf_setn(at,1);

    sf_oaxa(Fout,axcr,1);  
    sf_oaxa(Fout,at,2);  
    sf_oaxa(Fout,at,3);  

    //=========================== 

    // To adjoint source

    nr=sf_n(axcr); 

    /*-------------------------------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    /*-------------------------------------------------------------------------------------*/
    rr = pt2dalloc1(nr);
    pt2dread1(Fxcr,rr,nr,2);



    // allocate wavefield space    
    image = sf_floatalloc2(nz,nx);
    // read image 
    sf_floatread(image[0],nx*nz,Fimg);


    sparse_image = sf_floatalloc(nr);




    for (i1=0 ; i1 < nr ; i1++) {
    
        ix = 0.5+(rr[i1].x - ox)/dx;
        iz = 0.5+(rr[i1].z - oz)/dz;

        sparse_image[i1] = image[ix][iz] ;
    }

    sf_floatwrite(sparse_image,nr,Fout);

    
    free(*image); free(image);
    free(rr);
    free(sparse_image);
    exit(0);
}
