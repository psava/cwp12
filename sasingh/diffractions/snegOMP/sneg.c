/*  negentropy*/

#include <rsf.h>
#include <math.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
 
    int ompnth=1;


    int n1, n2,n3,box1,box2,klo1,khi1,klo2,khi2,kmid1,kmid2;
    int iy,i1,i2,d1,d2,o1,ix1,ix2,ix , iz,m1;
    float  ***d,***dat,***nsum,***neg,***dN;
    float sumG,neg1,AmpNorm,h,h1; 
    double scalea, logs;
    sf_file in, out; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in  = sf_input("in");
    /* standard output */
    out = sf_output("out");
    
#ifdef _OPENMP
   ompnth = omp_init();
#endif
    
    /* parameters from input file*/
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histint(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histint(in,"d2",&d2)) sf_error("No d2= in input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) n3=1;
fprintf(stderr,"values: %d %d\n",n1,n2,n3);
    /* parameter from the command line (i.e. box1=50 box2=50 ) */
    if (!sf_getint("box1",&box1)) sf_error("Need box1=");
    if (!sf_getint("box2",&box2)) sf_error("Need box2=");
    
    /* allocate floating point array */
    dat = sf_floatalloc3(n1,n2,n3);
    
    /* initialise the size of the searching box*/
    klo1=0;
    khi1=box1;
    klo2=0;
    khi2=box2;
    int s2= n2-(box2);
    int s1= n1-(box1);
    int bm1=box1/2;
    int bm2=box2/2;


    kmid1=((klo1+khi1)*0.5);
    kmid2=((klo2+khi2)*0.5)-1;
    nsum =sf_floatalloc3(n1,n2,n3);
    d    =sf_floatalloc3(n1,n2,n3);
    dN   =sf_floatalloc3(n1,n2,n3);
    neg  =sf_floatalloc3(n1,n2,n3);
    for (iy=0; iy<n3; ++iy){
      for (ix=0; ix<n2; ++ix) {
       for (iz=0; iz<n1; ++iz) {
          neg [iy][ix][iz]=0.0;
          nsum[iy][ix][iz]=0.0;
          dat [iy][ix][iz]=0.0;
          d   [iy][ix][iz]=0.0;
	        dN  [iy][ix][iz]=0.0;
        }
      }
     } 
    sf_floatread(dat[0][0],n1*n2*n3,in); 
    sumG=0;
        for(iy=0; iy<n3; ++iy){
          for (i2=o1; i2<n2; ++i2){
            for (i1=o1; i1<n1; ++i1){
	            d[iy][i2][i1]= dat[iy][i2][i1] * dat[iy][i2][i1]; //make all amplitudes positive
              sumG += d[iy][i2][i1];// sum of each trace
            }
	        }
        }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) \
    private(klo2,klo1,iy,kmid2,m1,khi1,khi2,neg1,AmpNorm,scalea,logs,h1,h) 
#endif

//initialise the box for negentropy
    for (iy=0; iy<n3; ++iy){  
      for (klo2=0; klo2<(s2) ; ++klo2) {
     
        for (klo1=0; klo1<(s1); ++klo1){
	
//  for (klo2=0, khi2=box2; klo2<(s2) && khi2<(n2); ++klo2,++khi2) {
//    kmid2++;
//    m1=kmid1-1;
//   for (klo1=0, khi1=box1; klo1<(s1) && khi1<(n1); ++klo1,++khi1){ 
//     m1++;
  
//#pragma omp barrier
    neg1=0.0;scalea=0.0;logs=0.0;AmpNorm=0.0;h=0,h1=0;
    khi2=klo2+box2; khi1=klo1+box1; kmid2=klo2+bm2; m1=klo1+bm1;  
          for (ix1=klo2; ix1<khi2; ++ix1){
            for (ix2=klo1; ix2<khi1; ++ix2){
          AmpNorm =(d[iy][ix1][ix2])/(sumG);
      
          dN[iy][ix1][ix2]=AmpNorm;
          h=(((ix2-m1)*(ix2-m1))+((ix1-kmid2)*(ix1-kmid2)))/(2*box1*box2*1.41); //penalty oper
          h=exp(-4*h);
          h1=(((box1*0.5)*(box1*0.5))+((box2*0.5)*(box2*0.5)))/(2*box1*box2*1.41); //penalty oper
          h1=exp(-4*h1);
            if (h1>= h) scalea=0;
              else
                scalea=box1*box2*AmpNorm*(h-h1);
            if (AmpNorm==0) {logs=0;}
//		else {logs=log(scalea);}
//		        logs=log(scalea);
		          else {	logs= scalea*scalea;}

			      neg1 += (scalea*logs);
		}
         }
        		neg[iy][kmid2][m1] = neg1/(box1*box2);
    }
    } 
     
    }

         sf_floatwrite(neg[0][0]  ,n1*n2*n3 ,out);  //write negentropy

    exit(0);
}
