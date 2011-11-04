/* Elastic time-domain FD modeling, automatically determines whether or not to use 3D or 2D, supports arbitrary types of anisotropy.
 *
 *
 * VERSION 1.0 - May 6, 2011
 *
 *
 * WARNING WARNING WARNING: **************
 * THIS CODE IMPLEMENTS THE ROTATED STAGGERED GRID DESCRIBED BY SAENGER(2000) 
 * THIS CODE IS NOT MEANT TO BE A PERMANENT REPLACEMENT,and is only a HOT FIX until a replacement code can be produced.
 *
 * Densities are now located 1/2 sample off the X and Z axes.  Additionally, you are to specify nx-1 and nz-1 densities,
 * where nx and nz are the number of CiJ points in the grid.
 *
 * The grid looks as follows:
 *
 *     o----o----o----o
 *     |              |
 *     | <>   <>   <> |
 *     |              |
 *     o----o----o----o
 *     |              |
 *     | <>   <>   <> |
 *     |              |
 *     o----o----o----o
 *     |              |
 *     | <>   <>   <> |
 *     |              |
 *     o----o----o----o
 *
 *     where o represents the locations where CIJ are specified.
 *     and <> represents locations where densities are specified.
 *     
 *     You should have: nx, ny, nz points for Cijs.
 *     You should have nx-1, ny-1, nz-1 points for densities.
 *RESUME DOCUMENTATION *********
 *

Elastic wave equation finite difference modeling in both 2D and 3D, using an explicit time-domain solver.

*** Please see the SConstruct in book/tutorial/ewe for a SConstruct that demonstrates how to use 
predefined functions for using this program. ***

This program is designed to be as generic as possible, and allows you to use files
with arbitrary models, and arbitrary source and receiver geometries.  Source types are
as generic as possible.  Supports arbitrary types of anisotropy as well.  

The downside to the generality, is that the program is not as performant as dedicated solvers
that are less flexible.  The program is parallelized using OpenMP, so be sure to use a compatible compiler to take
advantage of the performance boost.

===========  OPTIONS  ======================================= 
ani - The type of anisotropy for this simulation.  Valid options:
        For 2D:
        Orthorhombic = 0
        Triclinic = 1

        For 3D:
        Orthorhombic = 0
        Triclinic    = 1

        VTI, HTI, and Isotropic media are special cases of Orthorhombic media. 
        TTI media can be represented using Triclinic media.

cfl   - Execute the CFL check.  If the CFL check fails, then it will cause the program to fail. 
        The CFL check will check both the stability and accuracy conditions for both p-waves and
        s-waves. Depending on the type of anisotropy that you specify, the CFL condition will
        use a safety factor (that you can override if necessary).  
        
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

        For elastic FD, the P-wave most greatly influences the stability, as it moves fastest
        on the grid.  

        The stability condition gives us a LOWER bound on the grid sampling for a given dt.

        ------ ACCURACY -------
        The accuracy condition is related to the number of gridpoints per wavelength.  Thus,

        safety*vmin / fmax > N * sqrt(dx^2+dy^2+dz^2) 

        where vmin is the minimum wave velocity in the model (usually S-wave), fmax is some
        relative measure of the maximum frequency of your wavelet (usually 1.5*peak for Ricker), 
        N is the number of points desired per wavelength (5), and safety is a safety factor that 
        is dependent on the type of anisotropy.  

        For elastic FD, the S-wave most greatly impacts the accuracy of the solution, as the S-wave
        is typically much higher frequency and travels at slower wave speeds, meaning shorter 
        wavelengths.  

        The accuracy condition places an UPPER bound on the grid sampling.

        ---- SAFETY FACTOR -----
        The safety factor depends on the type of anisotropy specified, and attempts to place a lower
        bound on the slowest S-wave velocity (guess):

        Orthorhombic - (3/4)
        Triclinic    - (1/2)

        You can also override the safety factor using the safety parameter.
safety- Override the safety factor for the CFL condition.  This should be a floating point (0-1.0).
fmax  - An estimate of the highest frequency content in your wavelet (for Ricker use 1.5*peak)

fsrf  - Use a free surface at the top boundary (z=0).  
        WARNING: The free surface condition appears to introduce numerical artifacts into the simulation.  
        USE AT YOUR OWN RISK.

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
debug - Print debugging information.  This is more detailed than verbose.

srctype - An integer which determines where the source wavelet is injected
        in the simulation.  Valid options are:  
            0 - Acceleration source
            1 - Displacement source
            2 - Stress source
            3 - Tensor source
        The default option is 2: Acceleration source.
        For Stress, Displacement and Acceleration sources, your wavelet
        needs to have only 3 components (z,x,y).
        For a Tensor source, you must specify wavelet components for 
        all 3 (2D) or 6 (3D) tensor components in the following order:
        2D: tzz, txx, tzx
        3D: tzz, txx, tyy, tyz, tzx, txy

        Hint:  To inject an acoustic source, use a stress source,
            with equal components on all three components.

dabc  - Use a sponge layer to attenuate waves off the edge of the grid.  Use this in 
        combination with the nb parameter.
abcone- In addition to the sponge layer, using a severe ramp at the very edge of the expanded 
        sponge layer to severely attenuate zero-incidence waves at the boundaries. 
        It's not clear if this condition actually affects most computations.

nbell - Size of gaussian used to linearly interpolate curves.  A value of 5 seems to work well.  
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
nqz, nqx, oqz, oqx, nqy, oqy, - Allows you to set the parameters for the axes.  Leave as defaults.

=============BOUNDARY CONDITIONS ========================

This code enforces a fixed reflecting boundary condition at the 
edge of the computational domain.  The absorbing sponge is used
IN ADDITION to this condition.

=============FILE DESCRIPTIONS   ========================      

Fdat.rsf - An RSF file containing your data in the following format:
            axis 1 - source location
            axis 2 - wavefield component (z,x,y) order
            axis 3 - Time

Fwav.rsf - An RSF file containing your wavelet information.  For elastic modeling, the wavelet needs 
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
       NOTE: For tensor sources, you must have an appropriate number of components.  See srctype for more information.

cccc  - An N+1 dimensional RSF file that contains the values for the stiffness coefficients to be used 
           as the model for the modeling code.  So, for 2D, this would be a 3D array of values concatenated 
           together in the order as described in the anisotropy section.  Each coefficient file contains 
           the value of that coefficient for every point in space. 
           The axes for this file are: Axis 1: Z; Axis 2: X; Axis 3: Y;
    
        The stiffness tensor coefficients are defined uniformly as follows, where 
        --x---y---z--(y)-----(y) describes how the coefficients depend on space.
        |C11 C12 C13 C14 C15 C16|
        |    C22 C23 C24 C25 C26|
        |        C33 C34 C35 C36|
        |            C44 C45 C46|
        |                C55 C56|
        |                    C66|

        The tensor is assumed to be symmetric.  

        Order of the coefficients in the N+1 dimensional file...
        (First coefficient is the first 2D array in the 3D array).
        2D Anisotropy Modes:

        Orthorhombic: C11, C33, C55, C13
        "Triclinic:" C11, C13, C15, C33, C35, C55 
        ***Triclinic basically allows access to all coefs in 2D, but is not really triclinic media
        ------------------------------------------------------------
        (First coefficient is the first 3D array in the 4D array).
        3D Anisotropy Modes:

        Orthorhombic: C11, C22, C33, C44, C55, C66, C12, C13, C23
        Triclinic: C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, 
            C35, C36, C44, C45, C46, C55, C56, C66

   
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

wfl     - The name of the file to save the wavefield snapshots to.  This will be an N+2 
          dimensional file.  The file will be organized as follows:
              1-2(3) axes, spatial coordinates
              3(4) axis, wavefield components, in the Z,X,(Y) order
              4(5) axis, time, sequential snapshots
              ***The parentheses indicate what the axes will be for 3D models.

dat     - The name of the file to save the receiver data to.  The data has the format of:
	      spatial coordinates, then the data components of the elastic wavefield in the 
	      same order as the wavefield.  Lastly, time.

========== USEFUL COMMANDS  ============================= 	  

To view the wavefield snapshots (2D case):
sfwindow < Fwfl.rsf n3=1 f3=0 | sfgrey gainpanel=a pclip=100 | sfpen

To view the data (2D case):
sfwindow < Fdat.rsf n3=1 f3=0 | sfgrey gainpanel=a pclip=100 | sfpen

========== TROUBLESHOOTING ===============================

If you aren't getting output, or your output is full of Nans, make sure
that you have the proper dimensions for your wavelet files, and that
your input files make sense.

If your simulation passes the CFL check, but you are still getting
somewhat inaccurate (dispersive) results, then try changing the nbell
parameter.  Usually a value of nbell=5 works well, but sometimes you
may need to increase it to a value near nbell=10.  

Make sure your source and receiver points are located inside the 
model space, otherwise you will get all NaNs and the simulation will
run forever.

======= TIPS ========

If the simulation seems to slow down as it's running, its a pretty
good indication that the simulation has become unstable and is overflowing
with NaNs.


*/
/*
  Copyright (C) 2008 Colorado School of Mines
  
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

#include "fdutil.h"


#include <math.h>
#define NOP 3 /* derivative operator half-size */

/* Coefficients are:  (-1)^(k+1)/k * (n!)^2 / (n+k)!(n-k)!
for n = order of approximation i.e. eighth order
and k = coefficient number */
/* CENTERED derivatives 4th order */
#if 0
#define C1 +0.800000   /* +4/5    */
#define C2 -0.200000   /* -1/5    */
#define C3 +0.038095   /* +4/105  */
#define C4 -0.003571   /* -1/280  */
#endif
// These are the correct coefficients.
// We must use 2*normal coefficients because
// the derivatives are on half-grid spacing
#if 0
#define C1 +1.600000   /* +8/5    */
#define C2 -0.400000   /* -2/5    */
#define C3 +0.076190   /* +8/105  */
#define C4 -0.007142   /* -2/280  */
#endif
#if 0
// These work (8th-order stencils)
#define C1 +1.234091  /* +8/5    */
#define C2 -0.106650   /* -2/5    */
#define C3 +0.023036   /* +8/105  */
#define C4 -0.005342   /* -2/280  */
#define C5 +0.001077   /* -2/280  */
#define C6 -0.000166   /* -2/280  */
#define C7 +0.000017   /* -2/280  */
#define C8 -0.000001   /* -2/280  */
#endif 
#define C1 +1.171875  /* +8/5    */
#define C2 -0.065104   /* -2/5    */
#define C3 +0.004687   /* +8/105  */

#if 0
#define C1 +0.888889  /* +8/5    */
#define C2 -0.311111   /* -2/5    */
#define C3 +0.113131   /* +8/105  */
#define C4 -0.035354   /* -2/280  */
#define C5 +0.008702   /* -2/280  */
#define C6 -0.001554   /* -2/280  */
#define C7 +0.000178   /* -2/280  */
#define C8 -0.000010   /* -2/280  */
# endif

#if 0
#define C1 +1.33  /* +4/5    */
#define C2 -0.1666  /* -1/5    */
#endif
//#define C3 +0.076190   /* +4/105  */
//#define C4 -0.007142   /* -1/280  */



// Derivatives for strain to displacement grid
#if 0
#define DzH2(a,ix,iz,dz,dr) (1.0/dr)*(a[ix][iz]-a[ix-1][iz-1])
#define DxH2(a,ix,iz,dx,dr) (1.0/dr)*(a[ix][iz-1]-a[ix-1][iz])
#endif
#if 0
// 8 order
#define DzH2(a,ix,iz,dz,dr) (1.0/dr)*(       \
        C8*(a[ix+7][iz+7] - a[ix-8][iz-8]) + \
        C7*(a[ix+6][iz+6] - a[ix-7][iz-7]) + \
        C6*(a[ix+5][iz+5] - a[ix-6][iz-6]) + \
        C5*(a[ix+4][iz+4] - a[ix-5][iz-5]) + \
        C4*(a[ix+3][iz+3] - a[ix-4][iz-4]) + \
        C3*(a[ix+2][iz+2] - a[ix-3][iz-3]) + \
        C2*(a[ix+1][iz+1] - a[ix-2][iz-2]) + \
        C1*(a[ix][iz]     - a[ix-1][iz-1]) )
#define DxH2(a,ix,iz,dx,dr) (1.0/dr)*(       \
        C8*(a[ix+7][iz-8]-a[ix-8][iz+7])   + \
        C7*(a[ix+6][iz-7]-a[ix-7][iz+6])   + \
        C6*(a[ix+5][iz-6]-a[ix-6][iz+5])   + \
        C5*(a[ix+4][iz-5]-a[ix-5][iz+4])   + \
        C4*(a[ix+3][iz-4]-a[ix-4][iz+3])   + \
        C3*(a[ix+2][iz-3]-a[ix-3][iz+2])   + \
        C2*(a[ix+1][iz-2]-a[ix-2][iz+1])   + \
        C1*(a[ix][iz-1]-a[ix-1][iz])       )
#endif
#define DzH2(a,ix,iz,dz,dr) (1.0/dr)*(       \
        C3*(a[ix+2][iz+2] - a[ix-3][iz-3]) + \
        C2*(a[ix+1][iz+1] - a[ix-2][iz-2]) + \
        C1*(a[ix][iz]     - a[ix-1][iz-1]) )
#define DxH2(a,ix,iz,dx,dr) (1.0/dr)*(       \
        C3*(a[ix+2][iz-3]-a[ix-3][iz+2])   + \
        C2*(a[ix+1][iz-2]-a[ix-2][iz+1])   + \
        C1*(a[ix][iz-1]-a[ix-1][iz])       )

// These are fourth order derivatives, they work correctly but we 
// don't use them for consistency between 2D and 3D
#if 0
#define DzH2(a,ix,iz,dz,dr) (1.0/dr)*(       \
        C2*(a[ix+1][iz+1] - a[ix-2][iz-2]) + \
        C1*(a[ix][iz]     - a[ix-1][iz-1]) )
#define DxH2(a,ix,iz,dx,dr) (1.0/dr)*(       \
        C2*(a[ix+1][iz-2]-a[ix-2][iz+1])   + \
        C1*(a[ix][iz-1]-a[ix-1][iz])       )
#endif

#define Dz2(a,ix,iz,dx,dz,dr) (DzH2(a,ix,iz,dz,dr)-DxH2(a,ix,iz,dx,dr))*(dr/(2.0*dz))
#define Dx2(a,ix,iz,dx,dz,dr) (DzH2(a,ix,iz,dz,dr)+DxH2(a,ix,iz,dx,dr))*(dr/(2.0*dx))

// Derivatives for Displacement to Strain grid
#if 0 
#define DzH1(a,ix,iz,dz,dr) (1.0/dr)*(a[ix+1][iz+1]-a[ix][iz])
#define DxH1(a,ix,iz,dx,dr) (1.0/dr)*(a[ix+1][iz]-a[ix][iz+1])
#endif
//These dont work... Wrong coefficients
#if 0
#define DzH1(a,ix,iz,dz,dr) (1.0/dr)*(      \
        C8*(a[ix+8][iz+8]-a[ix-7][iz-7]) + \
        C7*(a[ix+7][iz+7]-a[ix-6][iz-6]) + \
        C6*(a[ix+6][iz+6]-a[ix-5][iz-5]) + \
        C5*(a[ix+5][iz+5]-a[ix-4][iz-4]) + \
        C4*(a[ix+4][iz+4]-a[ix-3][iz-3]) + \
        C3*(a[ix+3][iz+3]-a[ix-2][iz-2]) + \
        C2*(a[ix+2][iz+2]-a[ix-1][iz-1]) + \
        C1*(a[ix+1][iz+1]-a[ix][iz]) )
#define DxH1(a,ix,iz,dx,dr) (1.0/dr)*(   \
        C8*(a[ix+8][iz-7]-a[ix-7][iz+8]) + \
        C7*(a[ix+7][iz-6]-a[ix-6][iz+7]) + \
        C6*(a[ix+6][iz-5]-a[ix-5][iz+6]) + \
        C5*(a[ix+5][iz-4]-a[ix-4][iz+5]) + \
        C4*(a[ix+4][iz-3]-a[ix-3][iz+4]) + \
        C3*(a[ix+3][iz-2]-a[ix-2][iz+3]) + \
        C2*(a[ix+2][iz-1]-a[ix-1][iz+2]) + \
        C1*(a[ix+1][iz]-a[ix][iz+1])  )
#endif
#define DzH1(a,ix,iz,dz,dr) (1.0/dr)*(      \
        C3*(a[ix+3][iz+3]-a[ix-2][iz-2]) + \
        C2*(a[ix+2][iz+2]-a[ix-1][iz-1]) + \
        C1*(a[ix+1][iz+1]-a[ix][iz]) )
#define DxH1(a,ix,iz,dx,dr) (1.0/dr)*(   \
        C3*(a[ix+3][iz-2]-a[ix-2][iz+3]) + \
        C2*(a[ix+2][iz-1]-a[ix-1][iz+2]) + \
        C1*(a[ix+1][iz]-a[ix][iz+1])  )

// These are fourth order derivatives, they work correctly but we 
// don't use them for consistency between 2D and 3D
#if 0
#define DzH1(a,ix,iz,dz,dr) (1.0/dr)*(      \
        C2*(a[ix+2][iz+2]-a[ix-1][iz-1]) + \
        C1*(a[ix+1][iz+1]-a[ix][iz]) )
#define DxH1(a,ix,iz,dx,dr) (1.0/dr)*(   \
        C2*(a[ix+2][iz-1]-a[ix-1][iz+2]) + \
        C1*(a[ix+1][iz]-a[ix][iz+1])  )
#endif
#define Dz1(a,ix,iz,dx,dz,dr) (DzH1(a,ix,iz,dz,dr)-DxH1(a,ix,iz,dx,dr))*(dr/(2.0*dz))
#define Dx1(a,ix,iz,dx,dz,dr) (DzH1(a,ix,iz,dz,dr)+DxH1(a,ix,iz,dx,dr))*(dr/(2.0*dx))

#if 0
// old derviatives, these do not work correctly, do not use
#define Dx(a,ix,iz,s) (C4*(a[ix+4][iz] - a[ix-4][iz]) + \
		       C3*(a[ix+3][iz] - a[ix-3][iz]) +	  \
		       C2*(a[ix+2][iz] - a[ix-2][iz]) +		\
		       C1*(a[ix+1][iz] - a[ix-1][iz])  )*s
#define Dz(a,ix,iz,s) (C4*(a[ix][iz+4] - a[ix][iz-4]) + \
		       C3*(a[ix][iz+3] - a[ix][iz-3]) +	  \
		       C2*(a[ix][iz+2] - a[ix][iz-2]) +		\
		       C1*(a[ix][iz+1] - a[ix][iz-1])  )*s

#define Dx3(a,ix,iy,iz,s) (C4*(a[iy][ix+4][iz] - a[iy][ix-4][iz]) + \
			  C3*(a[iy][ix+3][iz] - a[iy][ix-3][iz]) +   \
			  C2*(a[iy][ix+2][iz] - a[iy][ix-2][iz]) +	\
			  C1*(a[iy][ix+1][iz] - a[iy][ix-1][iz])  )*s
#define Dy3(a,ix,iy,iz,s) (C4*(a[iy+4][ix][iz] - a[iy-4][ix][iz]) + \
			  C3*(a[iy+3][ix][iz] - a[iy-3][ix][iz]) +   \
			  C2*(a[iy+2][ix][iz] - a[iy-2][ix][iz]) +	\
			  C1*(a[iy+1][ix][iz] - a[iy-1][ix][iz])  )*s
#define Dz3(a,ix,iy,iz,s) (C4*(a[iy][ix][iz+4] - a[iy][ix][iz-4]) + \
			  C3*(a[iy][ix][iz+3] - a[iy][ix][iz-3]) +   \
			  C2*(a[iy][ix][iz+2] - a[iy][ix][iz-2]) +	\
			  C1*(a[iy][ix][iz+1] - a[iy][ix][iz-1])  )*s
#endif

// 3D derivatives.
//Derivatives for Displacement to Strain grid
// These are second order stencils that work

#if 0
// 8 order
#define D3_11(a,ix,iy,iz,dr) (1.0/dr)*( \
        C8*(a[iy+8][ix+8][iz+8]-a[iy-7][ix-7][iz-7]) + \
        C7*(a[iy+7][ix+7][iz+7]-a[iy-6][ix-6][iz-6]) + \
        C6*(a[iy+6][ix+6][iz+6]-a[iy-5][ix-5][iz-5]) + \
        C5*(a[iy+5][ix+5][iz+5]-a[iy-4][ix-4][iz-4]) + \
        C4*(a[iy+4][ix+4][iz+4]-a[iy-3][ix-3][iz-3]) + \
        C3*(a[iy+3][ix+3][iz+3]-a[iy-2][ix-2][iz-2]) + \
        C2*(a[iy+2][ix+2][iz+2]-a[iy-1][ix-1][iz-1]) + \
        C1*(a[iy+1][ix+1][iz+1]-a[iy][ix][iz]) )

#define D3_12(a,ix,iy,iz,dr) (1.0/dr)*(\
        C8*(a[iy+8][ix+8][iz-7]-a[iy-7][ix-7][iz+8]) + \
        C7*(a[iy+7][ix+7][iz-6]-a[iy-6][ix-6][iz+7]) + \
        C6*(a[iy+6][ix+6][iz-5]-a[iy-5][ix-5][iz+6]) + \
        C5*(a[iy+5][ix+5][iz-4]-a[iy-4][ix-4][iz+5]) + \
        C4*(a[iy+4][ix+4][iz-3]-a[iy-3][ix-3][iz+4]) + \
        C3*(a[iy+3][ix+3][iz-2]-a[iy-2][ix-2][iz+3]) + \
        C2*(a[iy+2][ix+2][iz-1]-a[iy-1][ix-1][iz+2]) + \
        C1*(a[iy+1][ix+1][iz]-a[iy][ix][iz+1]) )

#define D3_13(a,ix,iy,iz,dr) (1.0/dr)*(\
        C8*(a[iy-7][ix+8][iz+8]-a[iy+8][ix-7][iz-7]) + \
        C7*(a[iy-6][ix+7][iz+7]-a[iy+7][ix-6][iz-6]) + \
        C6*(a[iy-5][ix+6][iz+6]-a[iy+6][ix-5][iz-5]) + \
        C5*(a[iy-4][ix+5][iz+5]-a[iy+5][ix-4][iz-4]) + \
        C4*(a[iy-3][ix+4][iz+4]-a[iy+4][ix-3][iz-3]) + \
        C3*(a[iy-2][ix+3][iz+3]-a[iy+3][ix-2][iz-2]) + \
        C2*(a[iy-1][ix+2][iz+2]-a[iy+2][ix-1][iz-1]) + \
        C1*(a[iy][ix+1][iz+1]-a[iy+1][ix][iz]) )

#define D3_14(a,ix,iy,iz,dr) (1.0/dr)*(\
        C8*(a[iy-7][ix+8][iz-7]-a[iy+8][ix-7][iz+8]) + \
        C7*(a[iy-6][ix+7][iz-6]-a[iy+7][ix-6][iz+7]) + \
        C6*(a[iy-5][ix+6][iz-5]-a[iy+6][ix-5][iz+6]) + \
        C5*(a[iy-4][ix+5][iz-4]-a[iy+5][ix-4][iz+5]) + \
        C4*(a[iy-3][ix+4][iz-3]-a[iy+4][ix-3][iz+4]) + \
        C3*(a[iy-2][ix+3][iz-2]-a[iy+3][ix-2][iz+3]) + \
        C2*(a[iy-1][ix+2][iz-1]-a[iy+2][ix-1][iz+2]) + \
        C1*(a[iy][ix+1][iz]-a[iy+1][ix][iz+1]) )
#endif 
#define D3_11(a,ix,iy,iz,dr) (1.0/dr)*( \
        C3*(a[iy+3][ix+3][iz+3]-a[iy-2][ix-2][iz-2]) + \
        C2*(a[iy+2][ix+2][iz+2]-a[iy-1][ix-1][iz-1]) + \
        C1*(a[iy+1][ix+1][iz+1]-a[iy][ix][iz]) )

#define D3_12(a,ix,iy,iz,dr) (1.0/dr)*(\
        C3*(a[iy+3][ix+3][iz-2]-a[iy-2][ix-2][iz+3]) + \
        C2*(a[iy+2][ix+2][iz-1]-a[iy-1][ix-1][iz+2]) + \
        C1*(a[iy+1][ix+1][iz]-a[iy][ix][iz+1]) )

#define D3_13(a,ix,iy,iz,dr) (1.0/dr)*(\
        C3*(a[iy-2][ix+3][iz+3]-a[iy+3][ix-2][iz-2]) + \
        C2*(a[iy-1][ix+2][iz+2]-a[iy+2][ix-1][iz-1]) + \
        C1*(a[iy][ix+1][iz+1]-a[iy+1][ix][iz]) )

#define D3_14(a,ix,iy,iz,dr) (1.0/dr)*(\
        C3*(a[iy-2][ix+3][iz-2]-a[iy+3][ix-2][iz+3]) + \
        C2*(a[iy-1][ix+2][iz-1]-a[iy+2][ix-1][iz+2]) + \
        C1*(a[iy][ix+1][iz]-a[iy+1][ix][iz+1]) )
#if 0 
#define D3_11(a,ix,iy,iz,dr) (1.0/dr)*( \
        (a[iy+1][ix+1][iz+1]-a[iy][ix][iz]) )
#define D3_12(a,ix,iy,iz,dr) (1.0/dr)*(\
        (a[iy+1][ix+1][iz]-a[iy][ix][iz+1]) )
#define D3_13(a,ix,iy,iz,dr) (1.0/dr)*(\
        (a[iy][ix+1][iz+1]-a[iy+1][ix][iz]) )
#define D3_14(a,ix,iy,iz,dr) (1.0/dr)*(\
        (a[iy][ix+1][iz]-a[iy+1][ix][iz+1]) )
#endif

// These stencils work, but ARE SLOW.  
#define Dx3_1(a,ix,iy,iz,dx,dy,dz,dr) (dr/(4*dx))*(\
        D3_11(a,ix,iy,iz,dr)+D3_12(a,ix,iy,iz,dr)+D3_13(a,ix,iy,iz,dr)+D3_14(a,ix,iy,iz,dr) )
        
#define Dy3_1(a,ix,iy,iz,dx,dy,dz,dr) (dr/(4*dy))*(\
        D3_11(a,ix,iy,iz,dr)+D3_12(a,ix,iy,iz,dr)-D3_13(a,ix,iy,iz,dr)-D3_14(a,ix,iy,iz,dr) )

#define Dz3_1(a,ix,iy,iz,dx,dy,dz,dr) (dr/(4*dz))*(\
        D3_11(a,ix,iy,iz,dr)-D3_12(a,ix,iy,iz,dr)+D3_13(a,ix,iy,iz,dr)-D3_14(a,ix,iy,iz,dr) )
//Derivatives for Strain to Displacement grid
//
// These are second order derivatives that work
#if 0
#define D3_21(a,ix,iy,iz,dr) (1.0/dr)*( \
        C8*(a[iy+7][ix+7][iz+7]-a[iy-8][ix-8][iz-8]) + \
        C7*(a[iy+6][ix+6][iz+6]-a[iy-7][ix-7][iz-7]) + \
        C6*(a[iy+5][ix+5][iz+5]-a[iy-6][ix-6][iz-6]) + \
        C5*(a[iy+4][ix+4][iz+4]-a[iy-5][ix-5][iz-5]) + \
        C4*(a[iy+3][ix+3][iz+3]-a[iy-4][ix-4][iz-4]) + \
        C3*(a[iy+2][ix+2][iz+3]-a[iy-3][ix-3][iz-3]) + \
        C2*(a[iy+1][ix+1][iz+1]-a[iy-2][ix-2][iz-2]) + \
        C1*(a[iy][ix][iz]-a[iy-1][ix-1][iz-1]) )

#define D3_22(a,ix,iy,iz,dr) (1.0/dr)*(\
        C8*(a[iy+7][ix+7][iz-8]-a[iy-8][ix-8][iz+7]) + \
        C7*(a[iy+6][ix+6][iz-7]-a[iy-7][ix-7][iz+6]) + \
        C6*(a[iy+5][ix+5][iz-6]-a[iy-6][ix-6][iz+5]) + \
        C5*(a[iy+4][ix+4][iz-5]-a[iy-5][ix-5][iz+4]) + \
        C4*(a[iy+3][ix+3][iz-4]-a[iy-4][ix-4][iz+3]) + \
        C3*(a[iy+2][ix+2][iz-3]-a[iy-3][ix-3][iz+2]) + \
        C2*(a[iy+1][ix+1][iz-2]-a[iy-2][ix-2][iz+1]) + \
        C1*(a[iy][ix][iz-1]-a[iy-1][ix-1][iz]) )

#define D3_23(a,ix,iy,iz,dr) (1.0/dr)*(\
        C8*(a[iy-8][ix+7][iz+7]-a[iy+7][ix-8][iz-8]) + \
        C7*(a[iy-7][ix+6][iz+6]-a[iy+6][ix-7][iz-7]) + \
        C6*(a[iy-6][ix+5][iz+5]-a[iy+5][ix-6][iz-6]) + \
        C5*(a[iy-5][ix+4][iz+4]-a[iy+4][ix-5][iz-5]) + \
        C4*(a[iy-4][ix+3][iz+3]-a[iy+3][ix-4][iz-4]) + \
        C3*(a[iy-3][ix+2][iz+2]-a[iy+2][ix-3][iz-3]) + \
        C2*(a[iy-2][ix+1][iz+1]-a[iy+1][ix-2][iz-2]) + \
        C1*(a[iy-1][ix][iz]-a[iy][ix-1][iz-1]) )

#define D3_24(a,ix,iy,iz,dr) (1.0/dr)*(\
        C8*(a[iy-8][ix+7][iz-8]-a[iy+7][ix-8][iz+7])  + \
        C7*(a[iy-7][ix+6][iz-7]-a[iy+6][ix-7][iz+6])  + \
        C6*(a[iy-6][ix+5][iz-6]-a[iy+5][ix-6][iz+5])  + \
        C5*(a[iy-5][ix+4][iz-5]-a[iy+4][ix-5][iz+4])  + \
        C4*(a[iy-4][ix+3][iz-4]-a[iy+3][ix-4][iz+3])  + \
        C3*(a[iy-3][ix+2][iz-3]-a[iy+2][ix-3][iz+2])  + \
        C2*(a[iy-2][ix+1][iz-2]-a[iy+1][ix-2][iz+1])  + \
        C1*(a[iy-1][ix][iz-1]-a[iy][ix-1][iz]) )
#endif
#if 0
#define D3_21(a,ix,iy,iz,dr) (1.0/dr)*( \
        (a[iy][ix][iz]-a[iy-1][ix-1][iz-1]) )
#define D3_22(a,ix,iy,iz,dr) (1.0/dr)*(\
        (a[iy][ix][iz-1]-a[iy-1][ix-1][iz]) )
#define D3_23(a,ix,iy,iz,dr) (1.0/dr)*(\
        (a[iy-1][ix][iz]-a[iy][ix-1][iz-1]) )
#define D3_24(a,ix,iy,iz,dr) (1.0/dr)*(\
        (a[iy-1][ix][iz-1]-a[iy][ix-1][iz]) )
#endif
#define D3_21(a,ix,iy,iz,dr) (1.0/dr)*( \
        C3*(a[iy+2][ix+2][iz+3]-a[iy-3][ix-3][iz-3]) + \
        C2*(a[iy+1][ix+1][iz+1]-a[iy-2][ix-2][iz-2]) + \
        C1*(a[iy][ix][iz]-a[iy-1][ix-1][iz-1]) )

#define D3_22(a,ix,iy,iz,dr) (1.0/dr)*(\
        C3*(a[iy+2][ix+2][iz-3]-a[iy-3][ix-3][iz+2]) + \
        C2*(a[iy+1][ix+1][iz-2]-a[iy-2][ix-2][iz+1]) + \
        C1*(a[iy][ix][iz-1]-a[iy-1][ix-1][iz]) )

#define D3_23(a,ix,iy,iz,dr) (1.0/dr)*(\
        C3*(a[iy-3][ix+2][iz+2]-a[iy+2][ix-3][iz-3]) + \
        C2*(a[iy-2][ix+1][iz+1]-a[iy+1][ix-2][iz-2]) + \
        C1*(a[iy-1][ix][iz]-a[iy][ix-1][iz-1]) )

#define D3_24(a,ix,iy,iz,dr) (1.0/dr)*(\
        C3*(a[iy-3][ix+2][iz-3]-a[iy+2][ix-3][iz+2])  + \
        C2*(a[iy-2][ix+1][iz-2]-a[iy+1][ix-2][iz+1])  + \
        C1*(a[iy-1][ix][iz-1]-a[iy][ix-1][iz]) )

#define Dx3_2(a,ix,iy,iz,dx,dy,dz,dr) (dr/(4*dx))*(\
        D3_21(a,ix,iy,iz,dr)+D3_22(a,ix,iy,iz,dr)+D3_23(a,ix,iy,iz,dr)+D3_24(a,ix,iy,iz,dr) )
        
#define Dy3_2(a,ix,iy,iz,dx,dy,dz,dr) (dr/(4*dy))*(\
        D3_21(a,ix,iy,iz,dr)+D3_22(a,ix,iy,iz,dr)-D3_23(a,ix,iy,iz,dr)-D3_24(a,ix,iy,iz,dr) )

#define Dz3_2(a,ix,iy,iz,dx,dy,dz,dr) (dr/(4*dz))*(\
        D3_21(a,ix,iy,iz,dr)-D3_22(a,ix,iy,iz,dr)+D3_23(a,ix,iy,iz,dr)-D3_24(a,ix,iy,iz,dr) )

enum AniType {ORTHORHOMBIC=0,TRICLINIC=1};
enum SourceType {STRESS=2, DISPLACEMENT=1, ACCELERATION=0, TENSOR=3};

int main(int argc, char* argv[])
{
    /* Declare RSF params */
    bool verb,fsrf,snap,dabc,debug,abcone,cfl,is2d;
    int  jsnap,ntsnap,jdata,nbell;
    enum AniType type;
    enum SourceType srctype;
    float fmax,safety;
    
    int ompnth, ompchunk;

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fccc=NULL; /* velocity  */
    sf_file Fden=NULL; /* density   */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

    int     nt,nz,nx,ny=0,ns,nr,nc,nb;
    int     it,iz,ix,iy;
    float   dt,dz,dx,dy=0,idz,idx,idy;

    sf_axis at,ax,ay=NULL,az;
    sf_axis as,ar,ac;

    /* init RSF */
    sf_init(argc,argv);
   

   
    int tsrc;
    if(! sf_getint("srctype",&tsrc)) tsrc=0; /* source type, see comments */
    
    srctype = tsrc;

    /* MAKE SURE USER HAS DECLARED ANISOTROPY TYPE*/
    int ttype;
    if(! sf_getint("ani", &ttype)) ttype=-1;/* Anisotropy type, see comments */
    if(ttype == -1 || ttype >2 || ttype < 0) sf_error("Invalid Anisotropy type");
    type = ttype;

    /* I/O arrays */
    float***ww=NULL;           /* wavelet   */
    float **dd=NULL;           /* data      */
    /*------------------------------------------------------------*/
    /* Temp values for strain */
    float    szz,   sxx,   syy,   sxy,   syz,   szx;

    /* Temp values for cijs */
    float tc11, tc12, tc13, tc14, tc15, tc16;
    float       tc22, tc23, tc24, tc25, tc26;
    float             tc33, tc34, tc35, tc36;
    float                   tc44, tc45, tc46;
    float                         tc55, tc56;
    float                               tc66;
    /*------------------------------------------------------------*/
    /* execution flags */
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* use sponge absorbing BC */
    if(! sf_getbool("abcone",&abcone)) abcone=false; /* use sharp brake at end of boundary layer */
    if(! sf_getbool("debug",&debug)) debug=false; /* print debugging info */
    if(! sf_getbool("cfl",&cfl)) cfl=true; /* use CFL check, will cause program to fail if not satisfied */
    /*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("in" ); /* wavelet   */
    Fccc = sf_input ("ccc"); /* stiffness */
    Fden = sf_input ("den"); /* density   */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */

    /* print debugging information */
    /*------------------------------------------------------------*/
    /* Determine if 3D or 2D, test for existence of 4th axis, n=1 if none */
    sf_axis test = sf_iaxa(Fccc,4);
    if (sf_n(test) == 1) is2d = true;
    else is2d = false;
    
    
    if(!dabc){
        fprintf(stderr,"*********************\n");
        fprintf(stderr,"NOT USING ABSORBING BOUNDARY CONDITIONS!!!\n");
        fprintf(stderr,"THIS WILL CAUSE LARGE REFLECTIONS AT EDGES\n");
        fprintf(stderr,"TO TURN ON ABSORBING, SET DABC=y\n");
        fprintf(stderr,"*********************\n");
    }

    /* Make sure that we have the proper number of coefs in our stiffness file */
    sf_axis temp; int tempn;

    if (is2d) {
        temp = sf_iaxa(Fccc,3);
    } else {
        temp = sf_iaxa(Fccc,4);
    }
    tempn = sf_n(temp);
    sf_warning("Stiffness coefficient file has %d coefficients\n",tempn);

    switch(type) {
        case ORTHORHOMBIC:
            if (is2d && tempn != 4) sf_error("Orthorhombic requires 4 coefficients in 2D");
            else if ( !is2d && tempn != 9) sf_error("Orthorhombic requires 9 coefficients in 3D");
            break;
        case TRICLINIC:
            if (is2d && tempn != 6) sf_error("Triclinic requires 6 coefficients in 2D");
            else if ( !is2d && tempn != 21) sf_error("Triclinic requires 21 coefficients in 3D");
            break;
    }
    /*-------------------------------------------------------------*/
    /* Wavefield cut parameters */
    sf_axis   acz=NULL,acx=NULL,acy=NULL;
    int       nqz,nqx,nqy;
    float     oqz,oqx,oqy;
    float     dqz,dqx,dqy;

    /* axes */
    at = sf_iaxa(Fwav,3); sf_setlabel(at,"time"); 
    	sf_setunit(at,"s"); if(verb) sf_raxa(at); /* time */

    az = sf_iaxa(Fccc,1); sf_setlabel(az,"space z"); 
   	 sf_setunit(az,"km");if(verb) sf_raxa(az); /* depth */

    ax = sf_iaxa(Fccc,2); sf_setlabel(ax,"space x");  
    	sf_setunit(ax,"km");if(verb) sf_raxa(ax); /* space x */

    as = sf_iaxa(Fsou,2); sf_setlabel(as,"sources");  
    	sf_setunit(as,"km");if(verb) sf_raxa(as); /* sources */

    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"receivers");
    	sf_setunit(ar,"km");if(verb) sf_raxa(ar); /* receivers */

    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);

    if(!is2d){ 
    	ay = sf_iaxa(Fccc,3); sf_setlabel(ay,"space y");  
		sf_setunit(ay,"km");if(verb) sf_raxa(ay); /* space y */
    	ny = sf_n(ay); dy = sf_d(ay);
    }

#ifdef _OPENMP
    if (! sf_getint("ompnth",&ompnth)) ompnth = omp_get_num_threads();
    omp_set_num_threads(ompnth);
    if (! sf_getint("ompchunk",&ompchunk)) ompchunk=1;
    sf_warning("threads: %d chunk size: %d", ompnth,ompchunk);
#endif
   


    ns = sf_n(as);
    nr = sf_n(ar);

    if(snap){
        if(!sf_getint  ("nqz",&nqz)) nqz=sf_n(az);
        if(!sf_getint  ("nqx",&nqx)) nqx=sf_n(ax);

        if(!sf_getfloat("oqz",&oqz)) oqz=sf_o(az);
        if(!sf_getfloat("oqx",&oqx)) oqx=sf_o(ax);

        if(!sf_getfloat("dqz",&dqz)) dqz=sf_d(az);
        if(!sf_getfloat("dqx",&dqx)) dqx=sf_d(ax);

        if(!is2d){
            if(!sf_getint  ("nqy",&nqy)) nqy=sf_n(ay);
            if(!sf_getfloat("oqy",&oqy)) oqy=sf_o(ay);
        }
    }
   /*------------------------------------------------------------*/
    /* other execution parameters */
    if(! sf_getint("nbell",&nbell)) nbell=1;  /* bell size */
    if(verb) sf_warning("nbell=%d",nbell);
    if(! sf_getint("jdata",&jdata)) jdata=1;
    
    if( !sf_getint("nb",&nb)) nb=NOP; /* padding size for absorbing boundary */
//    if (nb < NOP) nb = NOP;

    
    if(snap) {  /* save wavefield every *jsnap* time steps */
	    if(! sf_getint("jsnap",&jsnap)) jsnap=nt;
    }
    /* CFL check parameters */
    if (cfl){
        if(! sf_getfloat("fmax",&fmax)) sf_error("CFL: you did not specify fmax");

        if (! sf_getfloat("safety",&safety) || safety < 0) {
            switch(type){
                case ORTHORHOMBIC:
                    safety = 0.75f;
                    break;
                case TRICLINIC:
                    safety = 0.50f;
                    break;
            }
        }
        sf_warning("CFL: safety margin %f",safety);
    }

    /* --------------2d code---------------- */
if (is2d){
    /* IO Arrays */
    pt2d   *ss=NULL;           /* sources   */
    pt2d   *rr=NULL;           /* receivers */
    /* FDM structure */
    fdm2d    fdm=NULL;
    abcone2d abcp=NULL,abcs=NULL;
    sponge   spo=NULL;

    /*------------------------------------------------------------*/
    float **tt=NULL;
    float **ro=NULL;           /* density   */
    float **c11=NULL;
    float **c33=NULL;
    float **c55=NULL;
    float **c13=NULL;
    float **c15=NULL;
    float **c35=NULL;
    float **vp,**vs;
    /*------------------------------------------------------------*/
    /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    float **umz,**uoz,**upz,**uaz,**utz; 
    float **umx,**uox,**upx,**uax,**utx;
    /* stress/strain tensor */ 
    float **tzz,**tzx,**txx; 
    /*------------------------------------------------------------*/
    /* linear interpolation weights/indices */
    lint2d cs,cr;
    /* wavefield cut params */
    float     **uc=NULL;
    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    fdm=fdutil_init(verb,fsrf,az,ax,nb,1);
    fdbell_init(nbell);

    // MODIFIED 
    // These coordinates define the grid that Ux and Uz are on.
    fdm->nzpad = sf_n(az)+2*nb+2*NOP-1;
    fdm->nxpad = sf_n(ax)+2*nb+2*NOP-1;
    fdm->nb    = nb+NOP;
    // These coordinates define the grid that Txx Tzz and Tzx are on.
    //
    fdm2d sfdm = fdutil_init(verb,fsrf,az,ax,nb,1);
    sfdm->nzpad = fdm->nzpad -1;
    sfdm->nxpad = fdm->nxpad -1;
    sfdm->nb    = nb+NOP;

    fdm2d dfdm = fdutil_init(verb,fsrf,az,ax,nb,1);
    dfdm->nx    = nx -1;
    dfdm->nz    = nz -1;
    dfdm->nzpad = fdm->nzpad;
    dfdm->nxpad = fdm->nxpad;
    dfdm->nb    =nb+NOP+1;

    sf_warning("new dimensions ux nz,nx = (%d,%d), txx nz,nx = (%d,%d)",
            fdm->nzpad,fdm->nxpad,sfdm->nzpad,sfdm->nxpad);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    /*------------------------------------------------------------*/
    /* 2D vector components */
    nc=2;
    ac=sf_maxa(nc,0,1);
    /*------------------------------------------------------------*/
    /* setup output data header */
    sf_oaxa(Fdat,ar,1);
    sf_oaxa(Fdat,ac,2);
    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,3);
    /* setup output wavefield header */
    if(snap) {

        acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
        acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
        /* TODO: check if the imaging window fits in the wavefield domain */

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
        sf_oaxa(Fwfl,ac, 3);
        sf_oaxa(Fwfl,at, 4);
    }

    if(debug) sf_warning("MADE OUTPUT WAVEFIELD INFO");

    /*------------------------------------------------------------*/
    /* source array */
    if(srctype == TENSOR) {
        ww=sf_floatalloc3(ns,3,nt); 
         sf_floatread(ww[0][0],nt*3*ns,Fwav);
    } else {
        ww=sf_floatalloc3(ns,nc,nt); 
        sf_floatread(ww[0][0],nt*nc*ns,Fwav);
    }

    /* data array */
    dd=sf_floatalloc2(nr,nc);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt2d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(nr,sizeof(*rr)); 

    pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
    pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

    cs = lint2d_make(ns,ss,fdm); /* setup linear interpolation */
    cr = lint2d_make(nr,rr,fdm);

    /*------------------------------------------------------------*/
    /* setup FD coefficients */

    idz = 1/dz;
    idx = 1/dx;

    float dr  = sqrt(dx*dx+dz*dz);


    if (debug) sf_warning("BEGINNING ALLOCATION");

    /*------------------------------------------------------------*/ 
    float **tro = sf_floatalloc2(nz-1,nx-1); 

    ro =sf_floatalloc2(fdm->nzpad,fdm->nxpad); /* input density */
    sf_floatread(tro[0],(nz-1)*(nx-1),Fden);     expand(tro,ro ,dfdm);
   
    free(*tro); free(tro);
    if (debug) sf_warning("ALLOCATED DENSITY");

    tt = sf_floatalloc2(nz,nx); 
    switch(type){
        case ORTHORHOMBIC:
            c11=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad);
            c33=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad);
            c55=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad);
            c13=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad);
            /* input stiffness */
            sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c11,sfdm);
            sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c33,sfdm);
            sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c55,sfdm);
            sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c13,sfdm);    
            break;
        case TRICLINIC:
            c11=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad); 
            c33=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad); 
            c55=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad); 
            c13=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad); 
            c15=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad);
            c35=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad);
            /* input stiffness */
            sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c11,sfdm);
            sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c13,sfdm); 
            sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c15,sfdm); 
            sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c33,sfdm);
            sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c35,sfdm);  
            sf_floatread(tt[0],nz*nx,Fccc );    expand(tt,c55,sfdm);              
            break;
    }
    free(*tt); free(tt);
    if (debug) sf_warning("DONE ALLOCATION");
    /*------------------------------------------------------------*/
	/* one-way abc setup   */
	vp = sf_floatalloc2(sfdm->nzpad,sfdm->nxpad); 
	vs = sf_floatalloc2(sfdm->nzpad,sfdm->nxpad); 
    float vpmax = 0.0f;
    float vpmin = 100000.0f;
    float vsmax = 0.0f;
    float vsmin = 100000.0f;
    float c11f;
    float c55f;
	for    (ix=0; ix<sfdm->nxpad; ix++) {
	    for(iz=0; iz<sfdm->nzpad; iz++) {
            c11f = c11[ix][iz]; 
            c55f = c55[ix][iz];
            if(c11f < 0 ) {
                vp[ix][iz] = sqrt( -c11f/ro[ix][iz] );
            } else {
                vp[ix][iz] = sqrt( c11f/ro[ix][iz] );
            }
            if (vp[ix][iz] > vpmax) vpmax = vp[ix][iz];
            if (vp[ix][iz] < vpmin) vpmin = vp[ix][iz];
            if(c55f < 0 ) {
                vs[ix][iz] = sqrt(-c55f/ro[ix][iz] );
            } else {
                vs[ix][iz] = sqrt( c55f/ro[ix][iz] );
            }
            if( vs[ix][iz] > vsmax) vsmax = vs[ix][iz];
            if( vs[ix][iz] < vsmin) vsmin = vs[ix][iz];
	    }
	}
    if (cfl) {
        cfl_elastic( vpmin,vpmax,vsmin,vsmax,
                     dx,-1.0,dz,dt,fmax,safety,4);
    }
    /* Absorbing boundary conditions setup */
    if (abcone) {
        abcp = abcone2d_make(NOP,dt,vp,fsrf,fdm);
        abcs = abcone2d_make(NOP,dt,vs,fsrf,fdm);
    }
	free(*vp); free(vp);
	free(*vs); free(vs);
	/* sponge abc setup */
    if (dabc){
	spo = sponge_make(fdm->nb);
    }

    /*------------------------------------------------------------*/
    /* precompute 1/ro * dt^2 */
    for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
            if (ro[ix][iz] < 1.0) sf_warning("ZERO DENSITY: ix %d iz %d",ix,iz);
	        ro[ix][iz] = dt*dt/ro[ix][iz];	    
	    }
    }
    /*------------------------------------------------------------*/
    if (debug) sf_warning("BEGIN ALLOCATION 2");
    /* allocate wavefield arrays */
    umz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uoz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    upz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uaz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    umx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uox=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    upx=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
    uax=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

    tzz=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad);
    tzx=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad);
    txx=sf_floatalloc2(sfdm->nzpad,sfdm->nxpad);

    /* Zero out values */
    if (debug) sf_warning("BEGIN ZEROING");

    for    (ix=0; ix<fdm->nxpad; ix++) {
	    for(iz=0; iz<fdm->nzpad; iz++) {
            umz[ix][iz]=0; umx[ix][iz]=0;
            uoz[ix][iz]=0; uox[ix][iz]=0;
            upz[ix][iz]=0; upx[ix][iz]=0;
            uaz[ix][iz]=0; uax[ix][iz]=0;
	    }
    }

    for    (ix=0; ix<sfdm->nxpad; ix++) {
	    for(iz=0; iz<sfdm->nzpad; iz++) {
            tzz[ix][iz]=0; tzx[ix][iz]=0; 
            txx[ix][iz]=0;
        }
    }

    if(debug) sf_warning("FINISHED ALL PRECOMPUTE");

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"%d/%d \r",it,nt);

	/*------------------------------------------------------------*/
	/* from displacement to strain                                */
	/*------------------------------------------------------------*/
	/* 
	 * exx = Fx(ux)
	 * ezz = Fz(uz)
	 * ezx = Bx(uz) + Bz(ux)
	 */
     if(debug) sf_warning("stepping stress");
#pragma omp parallel for	    \
    schedule(dynamic)		\
    private(iz,ix) \
    shared(sfdm,tzz,tzx,txx,uoz,uox,dx,dz,dr)
	for    (ix=NOP-1; ix < sfdm->nxpad-NOP+1; ix++) {
	    for(iz=NOP-1; iz < sfdm->nzpad-NOP+1; iz++) {
            txx[ix][iz] = Dx1(uox,ix,iz,dx,dz,dr);
	    }
	}		
#pragma omp parallel for	    \
    schedule(dynamic)		\
    private(iz,ix) \
    shared(sfdm,tzz,tzx,txx,uoz,uox,dx,dz,dr)
	for    (ix=NOP-1; ix < sfdm->nxpad-NOP+1; ix++) {
	    for(iz=NOP-1; iz < sfdm->nzpad-NOP+1; iz++) {
            tzz[ix][iz] = Dz1(uoz,ix,iz,dx,dz,dr);
        }
    }

#pragma omp parallel for	    \
    schedule(dynamic)		\
    private(iz,ix) \
    shared(sfdm,tzz,tzx,txx,uoz,uox,dx,dz,dr)
	for    (ix=NOP-1; ix < sfdm->nxpad-NOP+1; ix++) {
	    for(iz=NOP-1; iz < sfdm->nzpad-NOP+1; iz++) {
            tzx[ix][iz] = Dx1(uoz,ix,iz,dx,dz,dr) + Dz1(uox,ix,iz,dx,dz,dr);
        }
    }
	/*------------------------------------------------------------*/
	/* from strain to stress                                      */
	/*------------------------------------------------------------*/
    switch(type){
        case ORTHORHOMBIC:
#pragma omp parallel for	    \
    schedule(dynamic)		\
    private(iz,ix,szz,sxx,tc11,tc13,tc33)			\
    shared(sfdm,tzz,txx,c11,c33,c13)
            for    (ix=NOP-1; ix<sfdm->nxpad-NOP+1; ix++) {
                for(iz=NOP-1; iz<sfdm->nzpad-NOP+1; iz++) {
                    tc11 = c11[ix][iz];
                    tc13 = c13[ix][iz];
                    tc33 = c33[ix][iz];
                    
                    sxx = txx[ix][iz];
                    szz = tzz[ix][iz];

/*                    sxx = c11[ix][iz] * txx[ix][iz]
                        + c13[ix][iz] * tzz[ix][iz];
                            
                    szz = c13[ix][iz] * txx[ix][iz]
                        + c33[ix][iz] * tzz[ix][iz]; 
*/                    
                    txx[ix][iz] = tc11*sxx+tc13*szz;
                    tzz[ix][iz] = tc13*sxx+tc33*szz;
                }
            }
#pragma omp parallel for	    \
    schedule(dynamic)		\
    private(iz,ix,szx)			\
    shared(sfdm,tzx,c55)
            for    (ix=NOP-1; ix<sfdm->nxpad-NOP+1; ix++) {
                for(iz=NOP-1; iz<sfdm->nzpad-NOP+1; iz++) {
                    szx = c55[ix][iz] * tzx[ix][iz];
                    tzx[ix][iz] = szx;
                }
            }

        break;
        case TRICLINIC:
#pragma omp parallel for	    \
    schedule(dynamic)		\
    private(iz,ix,szz,szx,sxx)			\
    shared(sfdm,tzz,tzx,txx,c11,c33,c55,c13,c15,c35)
            for    (ix=NOP-1; ix<sfdm->nxpad-NOP+1; ix++) {
                for(iz=NOP-1; iz<sfdm->nzpad-NOP+1; iz++) {
                    sxx = c11[ix][iz] * txx[ix][iz]
                        + c13[ix][iz] * tzz[ix][iz]
                        + c15[ix][iz] * tzx[ix][iz];
                            
                    szz = c13[ix][iz] * txx[ix][iz]
                        + c33[ix][iz] * tzz[ix][iz] 
                        + c35[ix][iz] * tzx[ix][iz];
                    
                    szx = c15[ix][iz] * txx[ix][iz]
                        + c35[ix][iz] * tzz[ix][iz]
                        + c55[ix][iz] * tzx[ix][iz];

                    txx[ix][iz] = sxx;
                    tzz[ix][iz] = szz;

                    tzx[ix][iz] = szx;
                }
            }    
        break;
    }


	/*------------------------------------------------------------*/
	/* free surface */
	/*------------------------------------------------------------*/
	/* Francesco: the z component of the traction must be zero at the free surface */
	
    if(debug) sf_warning("freesurface");
	if(fsrf) {
        //sf_warning("USING BROKEN FREE SURFACE");
	    for(ix=sfdm->nb; ix < sfdm->nxpad-nb; ix++) {
		    for(iz=0; iz < sfdm->nb; iz++) {
//		        txx[ix][iz]=0;
		        tzz[ix][iz]=0;
		        tzx[ix][iz]=0;
		    }
/*            for(iz=sfdm->nb; iz < sfdm->nb+1; ++iz){
                tzz[ix][iz]=0;
                tzx[ix][iz]=0;
            }
            */
	    }
	}

	/*------------------------------------------------------------*/
	/* inject stress source                                       */
	/*------------------------------------------------------------*/
	if(srctype == STRESS || srctype == TENSOR) {
	    lint2d_bell(tzz,ww[it][0],cs);
	    lint2d_bell(txx,ww[it][1],cs);
	}

    if (srctype == TENSOR){
	    lint2d_bell(tzx,ww[it][2],cs);
    }
	
    if(debug) sf_warning("source");
    /*if(dabc){
        sponge2d_apply(txx,spo,fdm);
        sponge2d_apply(tzx,spo,fdm);
        sponge2d_apply(tzz,spo,fdm);    
	}*/

	/*------------------------------------------------------------*/
	/* from stress to acceleration                                */
	/*------------------------------------------------------------*/
	/* 
	 * ax = Bx(txx) + Fz(txz)
	 * az = Fx(txz) + Bz(tzz)
	 */
#pragma omp parallel for			\
    schedule(dynamic)		\
    private(iz,ix)				\
    shared(fdm,tzz,tzx,txx,uaz,uax,dx,dz,dr)
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
            uax[ix][iz] = Dx2( txx,ix,iz,dx,dz,dr) + Dz2( tzx,ix,iz,dx,dz,dr );
	    }
	}
 #pragma omp parallel for			\
    schedule(dynamic)		\
    private(iz,ix)				\
    shared(fdm,tzz,tzx,txx,uaz,uax,dx,dz,dr)
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
           uaz[ix][iz] = Dx2( tzx,ix,iz,dx,dz,dr) + Dz2( tzz,ix,iz,dx,dz,dr );
           }
    }

    if(debug) sf_warning("acceleration");
	/*------------------------------------------------------------*/
	/* inject acceleration source                                 */
	/*------------------------------------------------------------*/
	if(srctype == ACCELERATION) {
	    lint2d_bell(uaz,ww[it][0],cs); 
	    lint2d_bell(uax,ww[it][1],cs);
	}

	/*------------------------------------------------------------*/
	/* step forward in time                                       */
	/*------------------------------------------------------------*/
#pragma omp parallel for				\
    schedule(dynamic)			\
    private(iz,ix)					\
    shared(fdm,uoz,uox,umz,umx,upz,upx,uaz,uax,ro)
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
            upz[ix][iz] = 2*uoz[ix][iz] 
                -           umz[ix][iz] 
                +           uaz[ix][iz] * ro[ix][iz]; 

	    }
	}

#pragma omp parallel for				\
    schedule(dynamic) \
    private(iz,ix)					\
    shared(fdm,uoz,uox,umz,umx,upz,upx,uaz,uax,ro)
	for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
	    for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
            upx[ix][iz] = 2*uox[ix][iz] 
                -           umx[ix][iz] 
                +           uax[ix][iz] * ro[ix][iz]; 
	    }
	}

    if(debug) sf_warning("boundaries");
	if(srctype == DISPLACEMENT) {
	    lint2d_bell(upz,ww[it][0],cs); 
	    lint2d_bell(upx,ww[it][1],cs);
	}
	/*------------------------------------------------------------*/
	/* apply the boundary condition                               */
	/*------------------------------------------------------------*/
    /*
    for(ix = 0; ix < fdm->nxpad; ix++){
        for(iz = 0; iz < NOP; iz++){
            upz[ix][iz] = 0.0f;
            upx[ix][iz] = 0.0f;
        }
    }
    for(ix = 0; ix < fdm->nxpad; ix++){
        for(iz = fdm->nzpad-NOP; iz < fdm->nzpad; iz++){
            upz[ix][iz] = 0.0f;
            upx[ix][iz] = 0.0f;
        }
    }
    for(ix = 0; ix < NOP; ix++){
        for(iz = 0; iz < fdm->nzpad; iz++){
            upz[ix][iz] = 0.0f;
            upx[ix][iz] = 0.0f;
        }
    }
    for(ix = fdm->nxpad-NOP; ix < fdm->nxpad; ix++){
        for(iz = 0; iz < fdm->nzpad; iz++){
            upz[ix][iz] = 0.0f;
            upx[ix][iz] = 0.0f;
        }
    }
    */

    if(debug) sf_warning("sponging");
       /* sponge ABC */
    if(dabc){
        //sf_warning("I AM SPONGING EVERYWHERE");
        //sponge2d_apply(umz,spo,fdm);
        //sponge2d_apply(uoz,spo,fdm);
        sponge2d_apply(upz,spo,fdm);
        
        //sponge2d_apply(umx,spo,fdm);
        //sponge2d_apply(uox,spo,fdm);
        sponge2d_apply(upx,spo,fdm);
    }

	/* circulate wavefield arrays */
    /* Change pointers around */
	utz=umz; utx=umx;
	umz=uoz; umx=uox;
	uoz=upz; uox=upx;
	upz=utz; upx=utx;


    /* Apply the zero-incidence boundary condition */
    if(abcone){
            abcone2d_apply(uoz,umz,NOP,abcp,fdm);
            abcone2d_apply(uox,umx,NOP,abcp,fdm);
            
            abcone2d_apply(uoz,umz,NOP,abcs,fdm);
            abcone2d_apply(uox,umx,NOP,abcs,fdm);
    }
 	/*------------------------------------------------------------*/
	/* cut wavefield and save */
	/*------------------------------------------------------------*/
    if(debug) sf_warning("extracting");
	lint2d_extract(uoz,dd[0],cr);
	lint2d_extract(uox,dd[1],cr);

    /* Save snapshots of wavefield if selected */
	if(snap && it%jsnap==0) {
	    cut2d(uoz,uc,fdm,acz,acx);
	    sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
	    
	    cut2d(uox,uc,fdm,acz,acx);
	    sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
	}
    /* Save data snapshots*/
	if(it%jdata==0) sf_floatwrite(dd[0],nr*nc,Fdat);
    }

    if(verb) fprintf(stderr,"\n");    
    
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    
    free(**ww); free(*ww); free(ww);
    free(ss);
    free(rr);
    free(*dd);  free(dd);

    free(*ro);  free(ro);
    /* Only free double pointers that are not null... */
    switch (type){
        case TRICLINIC:
            free(*c15); free(c15);
            free(*c35); free(c35); //flow through
        case ORTHORHOMBIC:
            free(*c11); free(c11);
            free(*c33); free(c33);
            free(*c55); free(c55);
            free(*c13); free(c13);
        break;
    }
/* Free all outstanding parameters */
    free(*umz); free(umz);
    free(*uoz); free(uoz);
    free(*upz); free(upz);
    free(*uaz); free(uaz);

    free(*umx); free(umx);
    free(*uox); free(uox);
    free(*upx); free(upx);
    free(*uax); free(uax);

    free(*tzz); free(tzz);
    free(*txx); free(txx);
    free(*tzx); free(tzx);

    if(snap){free(*uc);  free(uc);    }
    exit (0);
/* ------------end 2d ------------------ */
} else {

/* ------------ 3d code ---------------*/  
    /* IO Arrays */
    pt3d   *ss=NULL;           /* sources   */
    pt3d   *rr=NULL;           /* receivers */
    /* FDM structure */
    fdm3d    fdm=NULL;
    abcone3d abcp=NULL,abcs=NULL;
    sponge   spo=NULL;
    /*------------------------------------------------------------*/
    float ***tt=NULL;
    float ***ro=NULL;           /* density */

    float ***c11=NULL,***c12=NULL,***c13=NULL,***c14=NULL,***c15=NULL,***c16=NULL;
    float ***c22=NULL,***c23=NULL,***c24=NULL,***c25=NULL,***c26=NULL;
    float ***c33=NULL,***c34=NULL,***c35=NULL,***c36=NULL;
    float ***c44=NULL,***c45=NULL,***c46=NULL;
    float ***c55=NULL,***c56=NULL;
    float ***c66=NULL;
  
    float ***vp,***vs;

    /*------------------------------------------------------------*/
    /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */
    float ***umz,***uoz,***upz,***uaz,***utz; 
    float ***umx,***uox,***upx,***uax,***utx;
    float ***umy,***uoy,***upy,***uay,***uty;
    /*------------------------------------------------------------*/

    /* stress/strain tensor */ 
    float ***tzz,***txx,***tyy,***txy,***tyz,***tzx;       

    /* linear interpolation weights/indices */
    lint3d cs,cr;

    /* wavefield cut params */
    float     ***uc=NULL;

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb)) nb=0;

    fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);
    fdbell3d_init(nbell);

    // MODIFIED 
    // These coordinates define the grid that Ux and Uz are on.
    fdm->nzpad = sf_n(az)+2*nb+2*NOP-1;
    fdm->nxpad = sf_n(ax)+2*nb+2*NOP-1;
    fdm->nypad = sf_n(ay)+2*nb+2*NOP-1;
    fdm->nb    = nb+NOP;
    // These coordinates define the grid that Txx Tzz and Tzx are on.
    //
    fdm3d sfdm = fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);
    sfdm->nzpad = fdm->nzpad -1;
    sfdm->nxpad = fdm->nxpad -1;
    sfdm->nypad = fdm->nypad -1;
    sfdm->nz = nz;
    sfdm->nx = nx;
    sfdm->ny = ny;
    sfdm->nb    = nb+NOP;
    sf_warning("new padded dimension sizes (%d,%d,%d) (z,x,y)",sfdm->nzpad,sfdm->nxpad,sfdm->nypad);
    sf_warning("new dimension sizes (%d,%d,%d) (z,x,y)",sfdm->nz,sfdm->nx,sfdm->ny);

    fdm3d dfdm = fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);
    dfdm->nx    = nx-1;
    dfdm->nz    = nz-1;
    dfdm->ny    = ny-1;
    dfdm->nzpad = fdm->nzpad;
    dfdm->nxpad = fdm->nxpad;
    dfdm->nypad = fdm->nypad;
    dfdm->nb    =nb+NOP+1;


    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); sf_setlabel(az,"expanded z");sf_setunit(az,"km");if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); sf_setlabel(ax,"expanded x");sf_setunit(ax,"km");if(verb) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); sf_setlabel(ay,"expanded y");sf_setunit(ay,"km");if(verb) sf_raxa(ay);
    /*------------------------------------------------------------*/

    /* 3D vector components */
    
    nc=3;
    ac=sf_maxa(nc,0,1);

    /*------------------------------------------------------------*/
    /* setup output data header */
    sf_oaxa(Fdat,ar,1);
    sf_oaxa(Fdat,ac,2);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,3);

    /* setup output wavefield header */
    if(debug) fprintf(stderr, "Setting up output wavefield header\n");
    if(snap) {
	dqz=sf_d(az);
	dqx=sf_d(ax);
	dqy=sf_d(ay);

	acz = sf_maxa(nqz,oqz,dqz);  sf_setunit(acz,"km");sf_setlabel(acz,"snapshot z");sf_raxa(acz);
	acx = sf_maxa(nqx,oqx,dqx);  sf_setunit(acx,"km");sf_setlabel(acx,"snapshot x");sf_raxa(acx);
	acy = sf_maxa(nqy,oqy,dqy);  sf_setunit(acy,"km");sf_setlabel(acy,"snapshot y");sf_raxa(acy);
	/* check if the imaging window fits in the wavefield domain */

	uc=sf_floatalloc3(sf_n(acz),sf_n(acx),sf_n(acy));

	ntsnap=0;
	for(it=0; it<nt; it++) {
	    if(it%jsnap==0) ntsnap++;
	}
	sf_setn(at,  ntsnap);
	sf_setd(at,dt*jsnap);
	sf_setlabel(at,"snapshot frames");
	if(verb)  sf_raxa(at);

	sf_oaxa(Fwfl,acz,1);
	sf_oaxa(Fwfl,acx,2);
	sf_oaxa(Fwfl,acy,3);
	sf_oaxa(Fwfl,ac, 4);
	sf_oaxa(Fwfl,at, 5);
    }

    /*------------------------------------------------------------*/
    /* source array */
    if (srctype == TENSOR) {
        ww = sf_floatalloc3(ns,6,nt);
        sf_floatread(ww[0][0],nt*6*ns,Fwav);
    }
    else {
        ww=sf_floatalloc3(ns,nc,nt); 
        sf_floatread(ww[0][0],nt*nc*ns,Fwav);
    }

    /* data array */
    dd=sf_floatalloc2(nr,nc);

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

    float dr = sqrt(dx*dx+dy*dy+dz*dz);


    /*------------------------------------------------------------*/ 
    tt = sf_floatalloc3(nz,nx,ny); 

    float ***tro = sf_floatalloc3(dfdm->nz,dfdm->nx,dfdm->ny);
    
    ro =sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    /* input density */
    sf_floatread(tro[0][0],(nz-1)*(nx-1)*(ny-1),Fden);     
    expand3d(tro,ro,dfdm);

    free(**tro); free(*tro); free(tro);

    if(debug) fprintf(stderr, "Beginning to read into stiffness arrays\n");

	/* allocate and read into arrays for stiffnesses */
    switch(type){
    	case ORTHORHOMBIC:
        if(debug)    sf_warning("allocating orthorhombic");
		c11=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
        if(debug)    sf_warning("allocating c11");
		c22=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c33=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c44=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c55=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c66=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c12=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c13=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c23=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad);  

        if(debug)    sf_warning("reading c11");
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c11,sfdm);
        if(debug)    sf_warning("expanding c11");
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c22,sfdm);    
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c33,sfdm);    
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c44,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c55,sfdm);    
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c66,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c12,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c13,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c23,sfdm);

	    break;

	case TRICLINIC:
		c11=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c12=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c13=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c14=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c15=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c16=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c22=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c23=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c24=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad);    
		c25=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c26=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c33=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c34=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c35=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c36=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c44=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c45=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c46=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad);    
		c55=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c56=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad); 
		c66=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad);  

		 /* input stiffness */
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c11,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c12,sfdm);    
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c13,sfdm);    
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c14,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c15,sfdm);    
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c16,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c22,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c23,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c24,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c25,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c26,sfdm);    
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c33,sfdm);    
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c34,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c35,sfdm);    
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c36,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c44,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c45,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c46,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c55,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c56,sfdm);
		sf_floatread(tt[0][0],nz*nx*ny,Fccc );    expand3d(tt,c66,sfdm);
	    break;
    }
    if(debug) fprintf(stderr,"Successfully allocated and read in coefficients...\n");
    free(**tt); free(*tt); free(tt);

    /*------------------------------------------------------------*/
    if(abcone) {
        /* one-way abc setup   */
        vp = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad); 
        vs = sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

        float vpmin = 100000000;
        float vpmax = 0.0f;
        float vsmin = 100000000;
        float vsmax = 0.0f;
        for (iy=0; iy<fdm->nypad; iy++) {
            for (ix=0; ix<fdm->nxpad; ix++) {
                for(iz=0; iz<fdm->nzpad; iz++) {
                    if (c11[iy][ix][iz] < 0) {
                        vp[iy][ix][iz] = sqrt(-c11[iy][ix][iz]/ro[iy][ix][iz]);
                    } else {
                        vp[iy][ix][iz] = sqrt(c11[iy][ix][iz]/ro[iy][ix][iz]);
                    } 
                    float vpt = vp[iy][ix][iz];
                    if (vpt > vpmax) vpmax = vpt;
                    if (vpt < vpmin) vpmin = vpt;
                    if (c55[iy][ix][iz] < 0){
                        vs[iy][ix][iz] = sqrt(-c55[iy][ix][iz]/ro[iy][ix][iz]);
                    } else {
                        vs[iy][ix][iz] = sqrt(c55[iy][ix][iz]/ro[iy][ix][iz]);
                    }
                    float vst = vs[iy][ix][iz];
                    if (vst > vsmax) vsmax = vst;
                    if (vst < vsmin) vsmin = vst;
                }
            }
        }

        if (cfl) {
            cfl_elastic(vpmin,vpmax,vsmin,vsmax,
                dx,dy,dz,dt,fmax,safety,4);
        }
        if (abcone){
            abcp = abcone3d_make(NOP,dt,vp,fsrf,fdm);
            abcs = abcone3d_make(NOP,dt,vs,fsrf,fdm);
        }
        free(**vp); free(*vp); free(vp);
        free(**vs); free(*vs); free(vs);
        /* sponge abc setup */
    }

    if(dabc){
        spo = sponge_make(fdm->nb);
    }
    /*------------------------------------------------------------*/
    /* precompute 1/ro * dt^2 */    
    for (iy=0; iy<fdm->nypad; iy++) {
	    for  (ix=0; ix<fdm->nxpad; ix++) {
	        for(iz=0; iz<fdm->nzpad; iz++) {
                if(ro[iy][ix][iz] < 1.0) sf_warning("DENSITY LESS THAN ONE (%d,%d,%d)",iz,ix,iy);
	    	    ro[iy][ix][iz] = dt*dt/ro[iy][ix][iz];
	        }
	    }
    }

    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
    if(debug) fprintf(stderr,"Allocating wavefield arrays\n");
    umz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uaz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    umx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uox=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upx=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uax=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    umy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uoy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    upy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
    uay=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);

    tzz=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad);
    tyy=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad);
    txx=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad);
    txy=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad);
    tyz=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad);
    tzx=sf_floatalloc3(sfdm->nzpad,sfdm->nxpad,sfdm->nypad);

    for (iy=0; iy<fdm->nypad; iy++) {
    	for (ix=0; ix<fdm->nxpad; ix++) {
	        for(iz=0; iz<fdm->nzpad; iz++) {
                umz[iy][ix][iz]=0; umx[iy][ix][iz]=0; umy[iy][ix][iz]=0;
                uoz[iy][ix][iz]=0; uox[iy][ix][iz]=0; uoy[iy][ix][iz]=0;
                upz[iy][ix][iz]=0; upx[iy][ix][iz]=0; upy[iy][ix][iz]=0;
                uaz[iy][ix][iz]=0; uax[iy][ix][iz]=0; uay[iy][ix][iz]=0;
	        }
	    }
    }

    for (iy=0; iy<sfdm->nypad; iy++) {
        for (ix=0; ix<sfdm->nxpad; ix++) {
	        for(iz=0; iz<sfdm->nzpad; iz++) {
                tzz[iy][ix][iz]=0; 
                tyy[iy][ix][iz]=0;
                txx[iy][ix][iz]=0;
                txy[iy][ix][iz]=0;
                tyz[iy][ix][iz]=0;
                tzx[iy][ix][iz]=0;
	        }
	    }
    }

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");

    for (it=0; it<nt; it++) {
	if(verb) fprintf(stderr,"%d/%d \r",it,nt);

	/*------------------------------------------------------------*/
	/* from displacement to strain                                */
	/*------------------------------------------------------------*/	
	/* 
	 * exx = Dx(ux)
	 * eyy = Dy(uy)
	 * ezz = Dz(uz)
	 * exy = Dy(ux) + Dx(uy)
	 * eyz = Dz(uy) + Dy(uz)
	 * ezx = Dx(uz) + Dz(ux)
	 */
#pragma omp parallel for					\
    schedule(dynamic)				\
    private(ix,iy,iz)						\
    shared(sfdm,txx,tyy,tzz,txy,tyz,tzx,uox,uoy,uoz,dx,dy,dz,dr)
	for        (iy=NOP-1; iy<sfdm->nypad-NOP+1; iy++) {
	    for    (ix=NOP-1; ix<sfdm->nxpad-NOP+1; ix++) {
    		for(iz=NOP-1; iz<sfdm->nzpad-NOP+1; iz++) { 
		    txx[iy][ix][iz] = Dx3_1(uox,ix,iy,iz,dx,dy,dz,dr);
            }
	    }
	}
#pragma omp parallel for					\
    schedule(dynamic)				\
    private(ix,iy,iz)						\
    shared(sfdm,txx,tyy,tzz,txy,tyz,tzx,uox,uoy,uoz,dx,dy,dz,dr)
	for        (iy=NOP-1; iy<sfdm->nypad-NOP+1; iy++) {
	    for    (ix=NOP-1; ix<sfdm->nxpad-NOP+1; ix++) {
    		for(iz=NOP-1; iz<sfdm->nzpad-NOP+1; iz++) { 
		    tyy[iy][ix][iz] = Dy3_1(uoy,ix,iy,iz,dx,dy,dz,dr);
            }
	    }
	}
#pragma omp parallel for					\
    schedule(dynamic)				\
    private(ix,iy,iz)						\
    shared(sfdm,txx,tyy,tzz,txy,tyz,tzx,uox,uoy,uoz,dx,dy,dz,dr)
	for        (iy=NOP-1; iy<sfdm->nypad-NOP+1; iy++) {
	    for    (ix=NOP-1; ix<sfdm->nxpad-NOP+1; ix++) {
    		for(iz=NOP-1; iz<sfdm->nzpad-NOP+1; iz++) { 
		    tzz[iy][ix][iz] = Dz3_1(uoz,ix,iy,iz,dx,dy,dz,dr);
            }
	    }
	}
#pragma omp parallel for					\
    schedule(dynamic)				\
    private(ix,iy,iz)						\
    shared(sfdm,txx,tyy,tzz,txy,tyz,tzx,uox,uoy,uoz,dx,dy,dz,dr)
	for        (iy=NOP-1; iy<sfdm->nypad-NOP+1; iy++) {
	    for    (ix=NOP-1; ix<sfdm->nxpad-NOP+1; ix++) {
    		for(iz=NOP-1; iz<sfdm->nzpad-NOP+1; iz++) { 
		    txy[iy][ix][iz] = Dy3_1(uox,ix,iy,iz,dx,dy,dz,dr) + Dx3_1(uoy,ix,iy,iz,dx,dy,dz,dr);
            }
	    }
	}
#pragma omp parallel for					\
    schedule(dynamic)				\
    private(ix,iy,iz)						\
    shared(sfdm,txx,tyy,tzz,txy,tyz,tzx,uox,uoy,uoz,dx,dy,dz,dr)
	for        (iy=NOP-1; iy<sfdm->nypad-NOP+1; iy++) {
	    for    (ix=NOP-1; ix<sfdm->nxpad-NOP+1; ix++) {
    		for(iz=NOP-1; iz<sfdm->nzpad-NOP+1; iz++) { 
		    tyz[iy][ix][iz] = Dz3_1(uoy,ix,iy,iz,dx,dy,dz,dr) + Dy3_1(uoz,ix,iy,iz,dx,dy,dz,dr);
            }
	    }
	}
#pragma omp parallel for					\
    schedule(dynamic)				\
    private(ix,iy,iz)						\
    shared(sfdm,txx,tyy,tzz,txy,tyz,tzx,uox,uoy,uoz,dx,dy,dz,dr)
	for        (iy=NOP-1; iy<sfdm->nypad-NOP+1; iy++) {
	    for    (ix=NOP-1; ix<sfdm->nxpad-NOP+1; ix++) {
    		for(iz=NOP-1; iz<sfdm->nzpad-NOP+1; iz++) { 
		    tzx[iy][ix][iz] = Dx3_1(uoz,ix,iy,iz,dx,dy,dz,dr) + Dz3_1(uox,ix,iy,iz,dx,dy,dz,dr);
            }
	    }
	}

	/*------------------------------------------------------------*/
	/* from strain to stress                                      */
	/*------------------------------------------------------------*/
if(debug) fprintf(stderr,"Going from strain to stress\n");
switch (type){

case TRICLINIC:
#pragma omp parallel for						\
    schedule(dynamic)					\
    private(ix,iy,iz,sxx,syy,szz,sxy,syz,szx)				\
    shared(sfdm,txx,tyy,tzz,txy,tyz,tzx,c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
	for        (iy=NOP-1; iy<sfdm->nypad-NOP+1; iy++) {
	    for    (ix=NOP-1; ix<sfdm->nxpad-NOP+1; ix++) {
		    for(iz=NOP-1; iz<sfdm->nzpad-NOP+1; iz++) {
		    
		    sxx = c11[iy][ix][iz] * txx[iy][ix][iz]
			+ c12[iy][ix][iz] * tyy[iy][ix][iz]
			+ c13[iy][ix][iz] * tzz[iy][ix][iz]
			+ c14[iy][ix][iz] * tyz[iy][ix][iz]
			+ c15[iy][ix][iz] * tzx[iy][ix][iz]
			+ c16[iy][ix][iz] * txy[iy][ix][iz];

		    syy = c12[iy][ix][iz] * txx[iy][ix][iz]
			+ c22[iy][ix][iz] * tyy[iy][ix][iz]
			+ c23[iy][ix][iz] * tzz[iy][ix][iz]
			+ c24[iy][ix][iz] * tyz[iy][ix][iz]
			+ c25[iy][ix][iz] * tzx[iy][ix][iz]
			+ c26[iy][ix][iz] * txy[iy][ix][iz];

		    szz = c13[iy][ix][iz] * txx[iy][ix][iz]
			+ c23[iy][ix][iz] * tyy[iy][ix][iz]
			+ c33[iy][ix][iz] * tzz[iy][ix][iz]
			+ c34[iy][ix][iz] * tyz[iy][ix][iz]
			+ c35[iy][ix][iz] * tzx[iy][ix][iz]
			+ c36[iy][ix][iz] * txy[iy][ix][iz];

		    syz = c14[iy][ix][iz] * txx[iy][ix][iz]
			+ c24[iy][ix][iz] * tyy[iy][ix][iz]
			+ c34[iy][ix][iz] * tzz[iy][ix][iz]
			+ c44[iy][ix][iz] * tyz[iy][ix][iz]
			+ c45[iy][ix][iz] * tzx[iy][ix][iz]
			+ c46[iy][ix][iz] * txy[iy][ix][iz];
		    
		    szx = c15[iy][ix][iz] * txx[iy][ix][iz]
			+ c25[iy][ix][iz] * tyy[iy][ix][iz]
			+ c35[iy][ix][iz] * tzz[iy][ix][iz]
			+ c45[iy][ix][iz] * tyz[iy][ix][iz]
			+ c55[iy][ix][iz] * tzx[iy][ix][iz]
			+ c56[iy][ix][iz] * txy[iy][ix][iz];

		    sxy = c16[iy][ix][iz] * txx[iy][ix][iz]
			+ c26[iy][ix][iz] * tyy[iy][ix][iz]
			+ c36[iy][ix][iz] * tzz[iy][ix][iz]
			+ c46[iy][ix][iz] * tyz[iy][ix][iz]
			+ c56[iy][ix][iz] * tzx[iy][ix][iz]
			+ c66[iy][ix][iz] * txy[iy][ix][iz];

		    txx[iy][ix][iz] = sxx;
		    tyy[iy][ix][iz] = syy;
		    tzz[iy][ix][iz] = szz;
		    txy[iy][ix][iz] = sxy;
		    tyz[iy][ix][iz] = syz;
		    tzx[iy][ix][iz] = szx;
		}
	    }
	}
break;

case ORTHORHOMBIC:
#pragma omp parallel for						\
    schedule(dynamic)					\
    private(ix,iy,iz,sxx,syy,szz,sxy,syz,szx)				\
    shared(sfdm,txx,tyy,tzz,txy,tyz,tzx,c11,c22,c33,c44,c55,c66,c12,c13,c23)
	for        (iy=NOP-1; iy<sfdm->nypad-NOP+1; iy++) {
	    for    (ix=NOP-1; ix<sfdm->nxpad-NOP+1; ix++) {
    		for(iz=NOP-1; iz<sfdm->nzpad-NOP+1; iz++) {
		    
		    sxx = c11[iy][ix][iz] * txx[iy][ix][iz]
			+ c12[iy][ix][iz] * tyy[iy][ix][iz]
			+ c13[iy][ix][iz] * tzz[iy][ix][iz];
		    syy = c12[iy][ix][iz] * txx[iy][ix][iz]
			+ c22[iy][ix][iz] * tyy[iy][ix][iz]
			+ c23[iy][ix][iz] * tzz[iy][ix][iz];
		    szz = c13[iy][ix][iz] * txx[iy][ix][iz]
			+ c23[iy][ix][iz] * tyy[iy][ix][iz]
			+ c33[iy][ix][iz] * tzz[iy][ix][iz];
		    
		    sxy = c66[iy][ix][iz] * txy[iy][ix][iz];
		    syz = c44[iy][ix][iz] * tyz[iy][ix][iz];
		    szx = c55[iy][ix][iz] * tzx[iy][ix][iz];
		    
		    txx[iy][ix][iz] = sxx;
		    tyy[iy][ix][iz] = syy;
		    tzz[iy][ix][iz] = szz;

		    txy[iy][ix][iz] = sxy;
		    tyz[iy][ix][iz] = syz;
		    tzx[iy][ix][iz] = szx;
		}
	    }
	}

break; }
    /* End Anisotropy LOGIC */
#if 0
   if(dabc) {
	    /* sponge ABC */
	    sponge3d_apply(txx,spo,fdm);
	    sponge3d_apply(tzx,spo,fdm);
	    sponge3d_apply(tyz,spo,fdm);
	    
	    sponge3d_apply(tzz,spo,fdm);
	    sponge3d_apply(tyy,spo,fdm);
	    sponge3d_apply(txy,spo,fdm);
	}
#endif
	/*------------------------------------------------------------*/
	/* free surface */
	/*------------------------------------------------------------*/
	/* Francesco: the z component of the traction must be zero at the free surface */
	if(fsrf) {
#if 0
#pragma omp parallel for						\
    schedule(static,fdm->ompchunk)					\
    private(ix,iy,iz)							\
    shared(sfdm,txx,tyy,tzz,txy,tyz,tzx)
#endif
	    for(iy=0; iy<sfdm->nypad; iy++) {
		    for(ix=0; ix<sfdm->nxpad; ix++) {
		        for(iz=0; iz < sfdm->nb; iz++) {
                    //txx[iy][ix][iz]=0;
                    //tyy[iy][ix][iz]=0;
                    tzz[iy][ix][iz]=0;
                    tyz[iy][ix][iz]=0;
                    tzx[iy][ix][iz]=0;
                    //txy[iy][ix][iz]=0;
                }
		    }
	    }
#if 0
#pragma omp parallel for						\
    schedule(static,fdm->ompchunk)					\
    private(ix,iy,iz)							\
    shared(sfdm,txx,tyy,tzz,txy,tyz,tzx)
#endif
	    for(iy=0; iy<sfdm->nypad; iy++) {
		    for(ix=0; ix<sfdm->nxpad; ix++) {
		        for(iz=fdm->nb; iz < sfdm->nb+1; iz++) {
                    tzz[iy][ix][iz]=0;
                    tyz[iy][ix][iz]=0;
                    tzx[iy][ix][iz]=0;
                }
            }
        }
	}

	/*------------------------------------------------------------*/
	/* inject stress source                                       */
	/*------------------------------------------------------------*/
	if(srctype == STRESS || srctype == TENSOR) {
	    lint3d_bell(txx,ww[it][0],cs);
	    lint3d_bell(tyy,ww[it][1],cs);
	    lint3d_bell(tzz,ww[it][2],cs);
	}
    if (srctype == TENSOR) {
        lint3d_bell(tyz,ww[it][3],cs);
        lint3d_bell(tzx,ww[it][4],cs);
        lint3d_bell(txy,ww[it][5],cs);
    }

	/*------------------------------------------------------------*/
	/* from stress to acceleration                                */
	/*------------------------------------------------------------*/
	/* 
	 * ax = Dx(txx) + Dy(txy) + Dz(txz)
	 * ay = Dx(txy) + Dy(tyy) + Dz(tyz)
	 * az = Dx(txz) + Dy(tyz) + Dz(tzz)
	 */	
     if(debug) fprintf(stderr,"Going from stress to acceleration \n");
#pragma omp parallel for					\
    schedule(dynamic)				\
    private(ix,iy,iz)						\
    shared(fdm,txx,tyy,tzz,txy,tyz,tzx,uax,uay,uaz,dx,dy,dz,dr)
	for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
		for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {		    
		    uax[iy][ix][iz] = Dx3_2( txx,ix,iy,iz,dx,dy,dz,dr) \
                            + Dy3_2( txy,ix,iy,iz,dx,dy,dz,dr) \
                            + Dz3_2( tzx,ix,iy,iz,dx,dy,dz,dr) ;
		    uay[iy][ix][iz] = Dx3_2( txy,ix,iy,iz,dx,dy,dz,dr) \
                            + Dy3_2( tyy,ix,iy,iz,dx,dy,dz,dr) \
                            + Dz3_2( tyz,ix,iy,iz,dx,dy,dz,dr) ;
		    uaz[iy][ix][iz] = Dx3_2( tzx,ix,iy,iz,dx,dy,dz,dr) \
                            + Dy3_2( tyz,ix,iy,iz,dx,dy,dz,dr) \
                            + Dz3_2( tzz,ix,iy,iz,dx,dy,dz,dr) ;		    
		}
	    }
	}

	/*------------------------------------------------------------*/
	/* inject acceleration source                                 */
	/*------------------------------------------------------------*/
	if(srctype == ACCELERATION) {
	    lint3d_bell(uaz,ww[it][0],cs);
	    lint3d_bell(uax,ww[it][1],cs);
	    lint3d_bell(uay,ww[it][2],cs);
	}

	/*------------------------------------------------------------*/
	/* step forward in time                                       */
	/*------------------------------------------------------------*/
    if(debug) fprintf(stderr,"Trying to step forward in time\n");
#pragma omp parallel for						\
    schedule(dynamic)					\
    private(ix,iy,iz)							\
    shared(fdm,uox,uoy,uoz,umx,umy,umz,upx,upy,upz,uax,uay,uaz,ro)
	for        (iy=NOP; iy<fdm->nypad-NOP; iy++) {
	    for    (ix=NOP; ix<fdm->nxpad-NOP; ix++) {
    		for(iz=NOP; iz<fdm->nzpad-NOP; iz++) {
                upx[iy][ix][iz] = 2*uox[iy][ix][iz] 
                -               umx[iy][ix][iz] 
                +               uax[iy][ix][iz] * ro[iy][ix][iz]; 

                upy[iy][ix][iz] = 2*uoy[iy][ix][iz] 
                -               umy[iy][ix][iz] 
                +               uay[iy][ix][iz] * ro[iy][ix][iz]; 

                upz[iy][ix][iz] = 2*uoz[iy][ix][iz] 
                -               umz[iy][ix][iz] 
                +               uaz[iy][ix][iz] * ro[iy][ix][iz]; 
                
            }
	    }
	}

	if(srctype == DISPLACEMENT) {
	    lint3d_bell(upz,ww[it][0],cs);
	    lint3d_bell(upx,ww[it][1],cs);
	    lint3d_bell(upy,ww[it][2],cs);
	}

	if(dabc) {
	    /* sponge ABC */
	    sponge3d_apply(umz,spo,fdm);
	    sponge3d_apply(uoz,spo,fdm);
	    sponge3d_apply(upz,spo,fdm);
	    
	    sponge3d_apply(umx,spo,fdm);
	    sponge3d_apply(uox,spo,fdm);
	    sponge3d_apply(upx,spo,fdm);

	    sponge3d_apply(umy,spo,fdm);
	    sponge3d_apply(uoy,spo,fdm);
	    sponge3d_apply(upy,spo,fdm);
	}	    


    if(debug) fprintf(stderr,"Done stepping forward\n");
	/* circulate wavefield arrays */
	utz=umz; uty=umy; utx=umx;
	umz=uoz; umy=uoy; umx=uox;
	uoz=upz; uoy=upy; uox=upx;
	upz=utz; upy=uty; upx=utx;
	
	    /* one-way ABC */
    if(abcone){
                abcone3d_apply(uoz,umz,NOP,abcp,fdm);
                abcone3d_apply(uox,umx,NOP,abcp,fdm);
                abcone3d_apply(uoy,umy,NOP,abcp,fdm);
                
                abcone3d_apply(uoz,umz,NOP,abcs,fdm);
                abcone3d_apply(uox,umx,NOP,abcs,fdm);
                abcone3d_apply(uoy,umy,NOP,abcs,fdm);
    }

	/*------------------------------------------------------------*/
	/* cut wavefield and save */
	/*------------------------------------------------------------*/
    if(debug) fprintf(stderr,"Trying to extract data\n");
	lint3d_extract(uoz,dd[0],cr);
	lint3d_extract(uox,dd[1],cr);
	lint3d_extract(uoy,dd[2],cr);
    if(debug) fprintf(stderr,"Trying to save snapshots\n");
	if(snap && it%jsnap==0) {
	    cut3d(uoz,uc,fdm,acz,acx,acy);
	    sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);

	    cut3d(uox,uc,fdm,acz,acx,acy);
	    sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);

	    cut3d(uoy,uc,fdm,acz,acx,acy);
	    sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
	}
	if(it%jdata==0) sf_floatwrite(dd[0],nr*nc,Fdat);
    }
    
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    if(debug) fprintf(stderr,"Finished loop, trying to deallocate...\n"); 
    free(**ww); free(*ww); free(ww);
    free(ss);
    free(rr);
    free(*dd);  free(dd);
    if(debug) fprintf(stderr,"Deallocating coefficient matrices...\n");

    free(**ro);  free(*ro);  free(ro);
    switch(type){
        case TRICLINIC:
            free(**c14); free(*c14); free(c14);
            free(**c15); free(*c15); free(c15);
            free(**c16); free(*c16); free(c16);
            free(**c24); free(*c24); free(c24);
            free(**c25); free(*c25); free(c25);
            free(**c26); free(*c26); free(c26);
            free(**c34); free(*c34); free(c34);
            free(**c35); free(*c35); free(c35);
            free(**c36); free(*c36); free(c36);
            free(**c45); free(*c45); free(c45);
            free(**c46); free(*c46); free(c46);
            free(**c56); free(*c56); free(c56);
        case ORTHORHOMBIC:
            free(**c11); free(*c11); free(c11);
            free(**c22); free(*c22); free(c22);
            free(**c33); free(*c33); free(c33);
            free(**c44); free(*c44); free(c44);
            free(**c55); free(*c55); free(c55);
            free(**c66); free(*c66); free(c66);
            free(**c12); free(*c12); free(c12);
            free(**c13); free(*c13); free(c13);
            free(**c23); free(*c23); free(c23);
            break;
    }

    free(**umz); free(*umz); free(umz);
    free(**uoz); free(*uoz); free(uoz);
    free(**upz); free(*upz); free(upz);
    free(**uaz); free(*uaz); free(uaz);

    free(**umx); free(*umx); free(umx);
    free(**uox); free(*uox); free(uox);
    free(**upx); free(*upx); free(upx);
    free(**uax); free(*uax); free(uax);

    free(**umy); free(*umy); free(umy);
    free(**uoy); free(*uoy); free(uoy);
    free(**upy); free(*upy); free(upy);
    free(**uay); free(*uay); free(uay);

    free(**tzz); free(*tzz); free(tzz);
    free(**txx); free(*txx); free(txx);
    free(**tyy); free(*tyy); free(tyy);
    free(**txy); free(*txy); free(txy);
    free(**tyz); free(*tyz); free(tyz);
    free(**tzx); free(*tzx); free(tzx);
    if(snap) {free(**uc);  free(*uc);  free(uc);}

    exit (0);
    }
}
/* ---------------- end 3d code --------------*/
