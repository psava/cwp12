/* Clip the data. */

#include <rsf.h>

int main(int argc, char* argv[])
{
    #define MIN(X,Y) ((X) < (Y) ? : (X) : (Y))
    #define MAX(X,Y) ((X) > (Y) ? : (X) : (Y))
    int n1, n2, i1, i2,i3, nlags, w_mean,n;
    int  k;
    float clip, *trace,*mean,*cov, scale;
    sf_file in, out, meanout; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in  = sf_input("in");
    /* standard output */
    out = sf_output("out");

    meanout = sf_output("mean");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype(in)) 
	sf_error("Need float input");

    /* n1 is the fastest dimension (trace length) */
    if (!sf_histint(in,"n1",&n1)) 
	sf_error("No n1= in input");
    /* leftsize gets n2*n3*n4*... (the number of traces) */
    n2 = sf_leftsize(in,1);

    /* half width of averaging operator*/
    if (!sf_getint("w_mean",&w_mean)) w_mean=10;

    /* lags for computting the covariance */
    if (!sf_getint("nlags",&nlags)) nlags=10; 

    if(nlags > n1) sf_error("lags cannot be greater that n1, exit now!");    

    scale= 1.0f/(2.0f*w_mean+1.0f);
    /* allocate floating point array */
    trace = sf_floatalloc (n1);
    mean  = sf_floatalloc (n1);
    cov   = sf_floatalloc (n1);


    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {

	/* read a trace */
	sf_floatread(trace,n1,in);
       
     
       /* I am assuming a mirror in the bounds */
 
	/* loop over samples to first mean */
	for (i1=-w_mean; i1 <  w_mean+1; i1++) {
             mean[0] += trace[abs(i1)];	
	}
        mean[0]*=scale;
        

	/* loop over samples to means without upper bound issues */
	for (i1=1; i1 < n1-w_mean; i1++) {
             mean[i1] = mean[i1-1] +scale*(trace[i1+w_mean] - trace[abs(i1-w_mean-1)]);	
	}

	/* loop over samples to means without upper bound issues */
	for (i1=n1-w_mean; i1 < n1; i1++) {
             i3= abs(n1-(i1+w_mean-n1)-1);
             mean[i1] = mean[i1-1] +scale*(trace[i3] - trace[abs(i1-w_mean-1)]);	
	}
	/* write a trace */
	sf_floatwrite(mean,n1,out);

        for (n=0;n<nlags; n++){
             cov[n]=0.0f;
             i1=0;   
             i3=i1+n;
             k=0;
             while (i3<n1){
                cov[n] += (trace[i1]-mean[i1])*(trace[i3]-mean[i3]);
		i1 +=1;
                i3 +=1;	
                k  +=1;
             }
             cov[n]*=1.0/k;
             fprintf(stderr,"%4d %10.4f\n",n,cov[n]);
         }
         sf_floatwrite(cov,nlags,meanout);

    }



    exit(0);
}
