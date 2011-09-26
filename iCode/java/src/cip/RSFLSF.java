//  RSFLSF.java
//
//  THOMAS CULLISON
//  Colorado School of Mines
//  Center for Wave Phenomena

package cip;

import java.util.*;
import static java.lang.Math.sqrt;
import java.io.IOException;
import edu.mines.jtk.dsp.LocalOrientFilter;
import edu.mines.jtk.dsp.LocalSemblanceFilter;

import rsf.RSF;
import rsf.Input;
import rsf.Output;

public class RSFLSF
{

	static { 
	  System.loadLibrary("jrsf");
	}

  enum SwitchOpt{
    LINSEM_2D,
    LINSEM_3D,
    PLNSEM_3D,
		DEFAULT
  };

	public static void main (String[] args){
		RSF parms = new RSF(args);
		Input input = new Input("in");
		Output output = new Output("out");

    String opt_str = parms.getString("opt","");
    if(opt_str.equals("")) {
			System.err.println("An option was not provided.");
			output.close();
			input.close();
			System.exit(-1);
		}

    int dim = parms.getInt("dim",0);
    if(dim <= 1 || 3 < dim) {
			System.err.println("The dimension must be equal to either 2 or 3");
			output.close();
			input.close();
			System.exit(-1);
		}

    SwitchOpt switch_opt = getSwitchOpt(opt_str,dim);

		//////////////////////////////////////////////////////
		// Get Semblance value zero-clip
		float zclip;
    zclip = parms.getFloat("zclip",0.f);
		assert(0.f <= zclip && zclip <= 1.f);

		//////////////////////////////////////////////////////
		// Get hw1 and hw2
		int hw1,hw2;
    hw1 = parms.getInt("hw1",0);
    if( hw1 <= 0 ) {
			System.err.println("The hw1 option must be greater than zero.");
			output.close();
			input.close();
			System.exit(-1);
		}
    hw2 = parms.getInt("hw2",0);
    if( hw2 < 0 ) {
			System.err.println("The hw2 option must be at least zero.");
			output.close();
			input.close();
			System.exit(-1);
		}


		//////////////////////////////////////////////////////
		// Get sigma(s)
		double sig1,sig2,sig3;
    sig1 = (double) parms.getFloat("sig1",0.f);
    if( sig1 <= 0.0 ) {
			System.err.println("The sig1 option must be greater than zero.");
			output.close();
			input.close();
			System.exit(-1);
		}
    sig2 = (double) parms.getFloat("sig2",0.f);
    if( sig2 <= 0.0 ) {
			System.err.println("The sig2 option has been set equal to sig1.");
			sig2 = sig1;
		}
		sig3 = (double) parms.getFloat("sig3",0.f);
		if( sig3 <= 0.0 ) {
			if(dim == 3) {
				System.err.println("The sig3 option has been set equal to sig1.");
				sig2 = sig1;
			}
		}


		//////////////////////////////////////////////////////
		// Get dims
		int n1,n2,n3;
		n1 = input.getN(1);
		n2 = input.getN(2);
		n3 = input.getN(3);


		//////////////////////////////////////////////////////////
		// Construct LocalOrientFilter and LocalSemblanceFilter
		LocalSemblanceFilter lsf = ImageLSFFactory.getLocalSemblanceFilter(hw1,hw2);
		LocalOrientFilter lof;
		if(dim == 3)
		  lof = ImageLOFFactory.getLocalOrientFilter(sig1,sig2,sig3);
		else
		  lof = ImageLOFFactory.getLocalOrientFilter(sig1,sig2);


		//////////////////////////////////////////////////////////
		// Carryout desired option
		float[][]   img2D;
		float[][][] img3D;
		float[][]   sem2D;
		float[][][] sem3D;
    switch (switch_opt)
    {
      case LINSEM_2D:
				img2D = new float[n2][n1];
				input.read(img2D);
				sem2D = ImageLSFFactory.getLinearSemblance2D(lof,lsf,img2D);
				if(zclip != 0.f)
				  zClip(sem2D,zclip);
				output.write(sem2D);
				output.close();
				input.close();
        break;
      case LINSEM_3D:
				img3D = new float[n3][n2][n1];
				input.read(img3D);
				sem3D = ImageLSFFactory.getLinearSemblance3D(lof,lsf,img3D);
				if(zclip != 0.f)
				  zClip(sem3D,zclip);
				output.write(sem3D);
				output.close();
				input.close();
        break;
      case PLNSEM_3D:
				img3D = new float[n3][n2][n1];
				input.read(img3D);
				sem3D = ImageLSFFactory.getPlanarSemblance3D(lof,lsf,img3D);
				if(zclip != 0.f)
				  zClip(sem3D,zclip);
				output.write(sem3D);
				output.close();
				input.close();
        break;
      default:
				System.err.println("Cannot find a match for the option provided");
				output.close();
				input.close();
				System.exit(-1);
    }

  } //end of main()


	/**********************************************************************/
	// 
	// Private Member Functions
	// 
	/**********************************************************************/

	//////////////////////////////////////////////////////
	// Translate argument option into an ENUM
  private static SwitchOpt getSwitchOpt(String sopt_str, int dim) {
		if(dim == 3) {
			if(sopt_str.equals("lin"))
				return SwitchOpt.LINSEM_3D;
			else if(sopt_str.equals("pln"))
				return SwitchOpt.PLNSEM_3D;
			else
				return SwitchOpt.DEFAULT;
		} 
		else {
			if(sopt_str.equals("lin"))
				return SwitchOpt.LINSEM_2D;
			else
				return SwitchOpt.DEFAULT;
		}
	}

	private static void zClip(float[][] sem, float zclip) {
	  int n1 = sem[0].length;
	  int n2 = sem.length;
		for(int i2=0; i2<n2; ++i2)
			for(int i1=0; i1<n1; ++i1)
			  if(sem[i2][i1] <= zclip)
					sem[i2][i1] = 0.f;
	}

	private static void zClip(float[][][] sem, float zclip) {
	  int n1 = sem[0][0].length;
	  int n2 = sem[0].length;
	  int n3 = sem.length;
		for(int i3=0; i3<n3; ++i3)
			for(int i2=0; i2<n2; ++i2)
				for(int i1=0; i1<n1; ++i1)
					if(sem[i3][i2][i1] <= zclip)
						sem[i3][i2][i1] = 0.f;
	}

}
