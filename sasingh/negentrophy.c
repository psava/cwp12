/*  negentropy*/

#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1, n2, box1,box2,klo1,khi1,klo2,khi2,kmid1,kmid2;
    float  **dat,**nsum,**neg,sumA;
    sf_file in, out; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in  = sf_input("in");
    /* standard output */
    out = sf_output("out");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype(in)) 
	sf_error("Need float input");

    /* n1 is the fastest dimension (trace length) */
    if (!sf_histint(in,"n1",&n1)) 
	sf_error("No n1= in input");
    /* leftsize gets n2*n3*n4*... (the number of traces) */
    n2 = sf_leftsize(in,1);

    /* parameter from the command line (i.e. box1=50 box2=50 ) */
    if (!sf_getfloat("box1",&box1)) sf_error("Need box1=");
    if (!sf_getfloat("box2",&box2)) sf_error("Nedd box2=");

    /* allocate floating point array */
    dat= sf_floatalloc2 (n1,n2);
	
    /* initialise the size of the searching box*/
    klo1=0;
    khi1=box1;
    klo2=0;
    khi2=box2;

    

    for (int klo1=0, int khi1=box1; khi1<=xmax; ++klo1,++khi1)
	
    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {

	/* read a trace */
	sf_floatread(trace,n1,in);

	/* loop over samples */
	for (i1=0; i1 < n1; i1++) {
	    if      (trace[i1] >  clip) trace[i1]= clip;
	    else if (trace[i1] < -clip) trace[i1]=-clip;
	}

	/* write a trace */
	sf_floatwrite(trace,n1,out);
    }


    exit(0);
}
