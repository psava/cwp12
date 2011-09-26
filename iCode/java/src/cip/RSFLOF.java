//  RSFLOF.java
//
//  THOMAS CULLISON
//  Colorado School of Mines
//  Center for Wave Phenomena

package cip;

import java.util.*;
import static java.lang.Math.sqrt;
import java.io.IOException;
import edu.mines.jtk.dsp.LocalOrientFilter;

//import rsf.Par;
import rsf.RSF;
import rsf.Input;
import rsf.Output;

public class RSFLOF
{

	static { 
	  System.loadLibrary("jrsf");
	}

  enum SwitchOpt{
    NORMAL_AT_COORD_2D,
    NORMAL_AT_COORD_3D,
    NORMAL_ALL_2D,
    NORMAL_ALL_3D,
    LINEARITY_2D,
    LINEARITY_3D,
		EIGENTENSOR_2D,
		EIGENTENSOR_3D,
    PLANARITY,
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
		// Get Linearity/Planarity value zero-clip
		float zclip;
    zclip = parms.getFloat("zclip",0.f);
		assert(0.f <= zclip && zclip <= 1.f);


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
				sig3 = sig1;
			}
		}


		//////////////////////////////////////////////////////
		// Get dims
		int n1,n2,n3;
		n1 = input.getN(1);
		n2 = input.getN(2);
		n3 = input.getN(3);


		//////////////////////////////////////////////////////
		// Construct LocalOrientFilter 
		LocalOrientFilter lof;
		if(dim == 3)
		  lof = ImageLOFFactory.getLocalOrientFilter(sig1,sig2,sig3);
		else
		  lof = ImageLOFFactory.getLocalOrientFilter(sig1,sig2);


		//////////////////////////////////////////////////////
		// Carryout desired option
		float[][]   img2D;
		float[][][] img3D;
		float[][]   lin2D;
		float[][][] lin3D;
		float[][][] pln3D;
		float[][][] normA2D;
		float[][][][] normA3D;
		float[][][] et2D;
		float[][][][] et3D;
		float[][] normC;
		float[][] coord;
		int npicks;
		float ox1,ox2,ox3,dx1,dx2,dx3;
		String cfile;
    switch (switch_opt)
    {
      case NORMAL_AT_COORD_2D:
				cfile = parms.getString("coordf","");
				if( cfile.equals("") ) {
					System.err.println("A coordinate file was not specified.");
					input.close();
					output.close();
					System.exit(-1);
				}
				Input cInput = new Input("coordf");
				npicks = cInput.getN(1);
				if( npicks == 1 )
					System.err.println("Warning!! Only one coordinate was specifed.");
				ox1 = 0.f;
				ox2 = 0.f;
				dx1 = 0.f;
				dx2 = 0.f;
				ox1 = parms.getFloat("o1",0.f);
				ox2 = parms.getFloat("o2",0.f);
				dx1 = parms.getFloat("d1",0.f);
				dx2 = parms.getFloat("d2",0.f);
				coord = new float[2][npicks];
				cInput.read(coord);
				normA2D = new float[2][n2][n1];
				input.read(normA2D);
				normC = ImageLOFFactory.getNormalByCoord2D(coord,normA2D,
																									 ox1, ox2, 
																									 dx1, dx2);
				output.setN(1,npicks);
				output.setN(2,2);
				output.setN(3,1);
				output.setN(4,1);
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);
				output.setOrigin(4,0.f);
				output.setDelta(1,0.f);
				output.setDelta(2,0.f);
				output.setDelta(3,0.f);
				output.setDelta(4,0.f);
				output.write(normC);
				input.close();
				cInput.close();
				output.close();
        break;
      case NORMAL_AT_COORD_3D:
				cfile = parms.getString("coordf","");
				if( cfile.equals("") ) {
					System.err.println("A coordinate file was not specified.");
					input.close();
					output.close();
					System.exit(-1);
				}
				Input cInput3D = new Input("coordf");
				npicks = cInput3D.getN(1);
				if( npicks == 1 )
					System.err.println("Warning!! Only one coordinate was specifed.");
				ox1 = 0.f;
				ox2 = 0.f;
				ox3 = 0.f;
				dx1 = 0.f;
				dx2 = 0.f;
				dx3 = 0.f;
				ox1 = parms.getFloat("o1",0.f);
				ox2 = parms.getFloat("o2",0.f);
				ox3 = parms.getFloat("o3",0.f);
				dx1 = parms.getFloat("d1",0.f);
				dx2 = parms.getFloat("d2",0.f);
				dx3 = parms.getFloat("d3",0.f);
				coord = new float[3][npicks];
				cInput3D.read(coord);
				normA3D = new float[3][n3][n2][n1];
				input.read(normA3D);
				normC = ImageLOFFactory.getNormalByCoord3D(coord, normA3D,
																									 ox1, ox2, ox3, 
																									 dx1, dx2, dx3);
				output.setN(1,npicks);
				output.setN(2,3);
				output.setN(3,1);
				output.setN(4,1);
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);
				output.setOrigin(4,0.f);
				output.setDelta(1,0.f);
				output.setDelta(2,0.f);
				output.setDelta(3,0.f);
				output.setDelta(4,0.f);
				output.write(normC);
				input.close();
				cInput3D.close();
				output.close();
        break;
      case NORMAL_ALL_2D:
				img2D = new float[n2][n1];
				input.read(img2D);
				normA2D = ImageLOFFactory.getNormalAll2D(lof,img2D);
				output.setN(3,2);
				output.setN(4,1);
				output.setOrigin(3,0.f);
				output.setOrigin(4,0.f);
				output.setDelta(3,0.f);
				output.setDelta(4,0.f);
				output.write(normA2D);
				output.close();
				input.close();
        break;
      case NORMAL_ALL_3D:
				img3D = new float[n3][n2][n1];
				input.read(img3D);
				normA3D = ImageLOFFactory.getNormalAll3D(lof,img3D);
				output.setN(4,3);
				output.setOrigin(4,0.f);
				output.setDelta(4,0.f);
				output.write(normA3D);
				output.close();
				input.close();
        break;
      case LINEARITY_2D:
				img2D = new float[n2][n1];
				input.read(img2D);
				lin2D = ImageLOFFactory.getLinearity2D(lof,img2D);
				if(zclip != 0.f)
				  zClip(lin2D,zclip);
				output.write(lin2D);
				output.close();
				input.close();
        break;
      case LINEARITY_3D:
				img3D = new float[n3][n2][n1];
				input.read(img3D);
				lin3D = ImageLOFFactory.getLinearity3D(lof,img3D);
				if(zclip != 0.f)
				  zClip(lin3D,zclip);
				output.write(lin3D);
				output.close();
				input.close();
        break;
      case PLANARITY:
				img3D = new float[n3][n2][n1];
				input.read(img3D);
				pln3D = ImageLOFFactory.getPlanarity3D(lof,img3D);
				if(zclip != 0.f)
				  zClip(pln3D,zclip);
				output.write(pln3D);
				output.close();
				input.close();
        break;
			/*
      case EIGENTENSOR_2D:
				img2D = new float[n2][n1];
				input.read(img2D);
				et2D = ImageLOFFactory.getEigenTensors2D(lof,img2D);
				output.setN(3,4);
				output.setOrigin(3,0.f);
				output.setDelta(3,0.f);
				output.write(et2D);
				output.close();
				input.close();
        break;
      case EIGENTENSOR_3D:
				img3D = new float[n3][n2][n1];
				input.read(img3D);
				et3D = ImageLOFFactory.getEigenTensors3D(lof,img3D);
				output.setN(4,9);
				output.setOrigin(4,0.f);
				output.setDelta(4,0.f);
				output.write(et3D);
				output.close();
				input.close();
        break;
			*/
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
				return SwitchOpt.LINEARITY_3D;
			if(sopt_str.equals("pln"))
				return SwitchOpt.PLANARITY;
			//if(sopt_str.equals("tensor"))
				//return SwitchOpt.EIGENTENSOR_3D;
			if(sopt_str.equals("nvec"))
				return SwitchOpt.NORMAL_ALL_3D;
			if(sopt_str.equals("natc"))
				return SwitchOpt.NORMAL_AT_COORD_3D;
			return SwitchOpt.DEFAULT;
		} 
		if(sopt_str.equals("lin"))
			return SwitchOpt.LINEARITY_2D;
		//if(sopt_str.equals("tensor"))
			//return SwitchOpt.EIGENTENSOR_2D;
		if(sopt_str.equals("nvec"))
			return SwitchOpt.NORMAL_ALL_2D;
		if(sopt_str.equals("natc"))
			return SwitchOpt.NORMAL_AT_COORD_2D;
		return SwitchOpt.DEFAULT;
	}

	private static void getOD(RSF parms,
														float ox1, float ox2,
	                          float dx1, float dx2) {

    ox1 = parms.getFloat("o1",0.f);
    ox2 = parms.getFloat("o2",0.f);
    dx1 = parms.getFloat("d1",0.f);
    dx2 = parms.getFloat("d2",0.f);
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


