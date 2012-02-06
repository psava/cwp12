/*  negentropy*/

#include <rsf.h>
#include <math.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
    int n1, n2, box1,box2,klo1,khi1,klo2,khi2,kmid1,kmid2;
    int i1,i2,d1,d2,o1,ix1,ix2,ix , iz,m1,ix3,ix4;

    float  **d,**dat,**nsum,**neg,**dN;
    float sumG,sumN,neg1,AmpNorm,h,h1; 
    double scalea, logs;
    sf_file in, out; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in  = sf_input("in");
    /* standard output */
    out = sf_output("out");
/*
    sf_file out2;
    out2= sf_output("data");
*/
    /* check that the input is float */
    //if (SF_FLOAT != sf_gettype(in)) 
	  //sf_error("Need float input");

    /* n1 is the fastest dimension (trace length) */
    //if (!sf_histint(in,"n1",&n1)) 
	  //sf_error("No n1= in input");
    /* leftsize gets n2*n3*n4*... (the number of traces) */
    //n2 = sf_leftsize(in,1);

    
    /* parameters from input file*/
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histint(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histint(in,"d2",&d2)) sf_error("No d2= in input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

    /* parameter from the command line (i.e. box1=50 box2=50 ) */
    if (!sf_getint("box1",&box1)) sf_error("Need box1=");
    if (!sf_getint("box2",&box2)) sf_error("Need box2=");
//    box1=box1*d1;
//    box2=box2*d2;
//    n1=n1*d1;
//    n2=n2*d2;
    fprintf(stderr,"values: %d %d %d %d\n",box1,box2,n1,n2);
    /* allocate floating point array */
    dat = sf_floatalloc2 (n1,n2);
	  //nsum= sf_floatalloc2 ((n1-
    /* initialise the size of the searching box*/
    klo1=0;
    khi1=box1;
    klo2=0;
    khi2=box2;
    int s1= n1-(box1);
    int s2= n2-(box2);
    fprintf(stderr,"values: %d %d\n",s1,s2);

    kmid1=((klo1+khi1)*0.5);
    kmid2=((klo2+khi2)*0.5)-1;

    fprintf(stderr,"values: %d %d\n",kmid1,kmid2);
    fprintf(stderr,"values: %d %d\n",n1,n2);

    fprintf(stderr,"I GOT HERE\n");
    nsum =sf_floatalloc2 (n1,n2);
    d    =sf_floatalloc2 (n1,n2);
    dN   =sf_floatalloc2 (n1,n2);
    neg  =sf_floatalloc2 (n1,n2);
    for (ix=0; ix<n2; ix++) {
      for (iz=0; iz<n1; iz++) {
        neg[ix][iz]=0.0;
        nsum[ix][iz]=0.0;
        dat[ix][iz]=0.0;
        d[ix][iz]=0.0;
	      dN[ix][iz]=0.0;
      }
    }
    fprintf(stderr,"I GOT HERE\n");
    sf_floatread(dat[0],n1*n2,in); 
//    sf_floatwrite(dat[0],n1*n2,out);
sumG=0;
    fprintf(stderr,"values: %d %d\n",n1,n2);
          for (i2=o1; i2<n2; ++i2){
            for (i1=o1; i1<n1; ++i1){
//              fprintf(stderr, "values: %f \n",dat[ix1][ix2]);
                if (dat[i2][i1] * dat[i2][i1] > (0.5e-20)) {
	            d[i2][i1]= dat[i2][i1] * dat[i2][i1]; //make all amplitudes positive
//              fprintf(stderr, "values1: %f \n",d[ix1][ix2]);
                  }
                  else d[i2][i1]=0;
//    fprintf(stderr,"values: %d %d\n",kmid1,kmid2);
              sumG += d[i2][i1];// sum of each trace
            }
	}
	fprintf(stderr,"values1: %f %f\n",d[50][50],sumG);
//    fprintf(stderr,"values: %d %d %d %d %f\n",ix1,m1,ix2,kmid1,sumA);
    //initialise the box for negentropy
    fprintf(stderr,"I GOT HERE1\n");
   for (klo2=0, khi2=box2; klo2<(n2-box2) && khi2<(n2); ++klo2,++khi2) {
        kmid2++; 
      m1=(kmid1-1);
//      fprintf(stderr,"values: %d %d \n", kmid1, m1); 
      for (klo1=0, khi1=box1; klo1<(n1-box1) && khi1<(n1); ++klo1,++khi1){
        m1++;
//    fprintf(stderr,"I GOT HERE LOOP\n");
//    fprintf(stderr,"values: %d %d %d %d %d\n",khi1,khi2,m,n1);
//          sumA=0;
	          sumN=0;
          for (ix3=klo2; ix3<khi2; ++ix3){
            for (ix4=klo1; ix4<khi1; ++ix4){
              sumN += d[ix3][ix4];}
              }
              
		neg1=0.0;scalea=0.0;logs=0.0;AmpNorm=0.0;h=0;
          for (ix1=klo2; ix1<khi2; ++ix1){
            for (ix2=klo1; ix2<khi1; ++ix2){
//              fprintf(stderr, "values: %f \n",dat[ix1][ix2]);
//	            d[ix1][ix2]= dat[ix1][ix2] * dat[ix1][ix2]; //make all amplitudes positive
//              fprintf(stderr, "values1: %f \n",d[ix1][ix2]);

//    fprintf(stderr,"values: %d %d\n",kmid1,kmid2);
      if (sumN==0.) AmpNorm=0;
      else (AmpNorm =(d[ix1][ix2])/(SumN));
       
	      dN[ix1][ix2]=AmpNorm;
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
//    fprintf(stderr,"values: %d %d %f %f\n",m1,ix2,sumA,nsum[m1][ix2]);
         }
		neg[kmid2][m1] = neg1/(box1*box2);
//    fprintf(stderr,"I GOT HERE LOOP\n");
         //normalising the amplitudes
//          logs= (scalea)*scalea;
	
//          logs= scalea;
//          neg[kmid1][kmid2] += scalea/box2;
//    fprintf(stderr,"values: %d %d %d\n",kmid1,m1,kmid2);
//          fprintf(stderr,"values: %f\n",neg);
//    fprintf(stderr,"values: %d %d %d\n",kmid1,m1,kmid2);

//         sf_floatwrite(neg[kmid1][kmid2],s1*s2,out);  //write negentropy
//         sf_floatwrite(neg[0],n1*n2,out);  //write negentropy
//         sf_floatwrite(nsum[0],n1*n2,out);  //write negentropy
//    fprintf(stderr,"values: %d %d\n",kmid1,kmid2);
    }
    
//    fprintf(stderr,"values: %d %d\n",kmid1,kmid2);
    }

//         sf_floatwrite(neg[0],n1*n2 ,out);  //write negentropy
         sf_floatwrite(neg[0]  ,n1*n2 ,out);  //write negentropy
//    fprintf(stderr,"I GOT HERE\n");
//    sf_floatwrite(neg,

    exit(0);
}
