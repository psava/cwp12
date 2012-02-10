/* adds all elements in a imagei*/


#include <rsf.h>
#include <math.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
    int n1, n2,o1,d1,d2,i1,i2,ix,iz; 

    float  ***dat,**sum;
    float sum1; 
    sf_file in, out; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in  = sf_input("in");
    /* standard output */
    out = sf_output("out");
    
    /* parameters from input file*/
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histint(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histint(in,"d2",&d2)) sf_error("No d2= in input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    fprintf(stderr,"values: %d %d\n",n1,n2,n3);
   sum = sf_floatalloc2(n3,n3);
     dat = sf_floatalloc3 (n1,n2,n3);

    fprintf(stderr,"values: %d %d\n",n1,n2);
  dat = sf_floatalloc2 (n1,n2);

   for (ix=0; ix<n2; ++ix) {
      for (iz=0; iz<n1; ++iz) {

           dat[ix][iz]=0.0;
      }
     }
sum1=0;

    fprintf(stderr,"values1: %f\n",sum1);
    sf_floatread(dat[0],n1*n2,in); 

//sum1=0;

    fprintf(stderr,"values1: %f\n",sum1);
      for (i2=o1; i2<n2; ++i2){
            for (i1=o1; i1<n1; ++i1){
      sum1 += dat[i2][i1];
        }
        }
    fprintf(stderr,"values: %f\n",sum1);
   sum[0][0]=sum1; 
    sf_floatwrite(sum[0],1*1,out);
      

  exit(0);
}



