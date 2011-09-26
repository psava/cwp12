/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package cip;

import java.util.*;
import static java.lang.Math.sqrt;
import java.io.IOException;
import edu.mines.jtk.dsp.LocalOrientFilter;
import edu.mines.jtk.dsp.LocalSemblanceFilter;

//import rsf.Par;
import rsf.RSF;
import rsf.Input;
import rsf.Output;

/**
 * Creates and writes EigenTensor(2|3) files in RSF file format.
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.11.14
 */
public class RSFTensor
{

	static { 
	  System.loadLibrary("jrsf");
	}

  enum SwitchOpt{
    METRIC_2D,
    METRIC_3D,
    STRUCTURE_2D,
    STRUCTURE_3D,
    ORIENTED_2D,
    ORIENTED_3D,
    EXCLUSION_2D,
    EXCLUSION_3D,
		TENSOR_AT_COORD_2D,
		TENSOR_AT_COORD_3D,
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
		int hw1,hw2;
		float ox1,ox2,ox3,dx1,dx2,dx3;
		String cfile;
		LocalSemblanceFilter lsf;

    switch (switch_opt)
    {
      case METRIC_2D:
				System.err.println("RSFTensor: 2D METRIC SEMBLANCE TENSORS");
				//////////////////////////////////////////////////////
				// Get hw1 and hw2
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
				lsf = new LocalSemblanceFilter(hw1,hw2);

				/////////////////////////////////////////////////////////////////
				// Read image file
				img2D = new float[n2][n1];
				System.err.println("Reading image file.");
				input.read(img2D);

				/////////////////////////////////////////////////////////////////
				// Compute tensors
				System.err.println("Computing metric-semblance tensors.");
				et2D = ImageTensorFactory.getMetricTensors2D(lof,lsf,img2D);

				/////////////////////////////////////////////////////////////////
				// Write tensors file and close inpute
				System.err.println("Writing metric-semblance-tensors file.");
				output.setN(3,4);
				output.setOrigin(3,0.f);
				output.setDelta(3,0.f);
				output.write(et2D);
				output.close();
				input.close();
        break;
      case METRIC_3D:
				System.err.println("RSFTensor: 3D METRIC SEMBLANCE TENSORS");
				//////////////////////////////////////////////////////
				// Get hw1 and hw2
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
				lsf = new LocalSemblanceFilter(hw1,hw2);

				/////////////////////////////////////////////////////////////////
				// Read image file
				img3D = new float[n3][n2][n1];
				System.err.println("Reading image file.");
				input.read(img3D);

				/////////////////////////////////////////////////////////////////
				// Compute tensors
				System.err.println("Computing metric-semblance tensors.");
				et3D = ImageTensorFactory.getMetricTensors3D(lof,lsf,img3D);

				/////////////////////////////////////////////////////////////////
				// Write tensors file and close inpute
				System.err.println("Writing metric-semblance-tensors file.");
				output.setN(4,9);
				output.setOrigin(4,0.f);
				output.setDelta(4,0.f);
				output.write(et3D);
				output.close();
				input.close();
        break;
      case STRUCTURE_2D:
				System.err.println("RSFTensor: 2D STRUCTURE TENSORS");

				/////////////////////////////////////////////////////////////////
				// Read image file
				img2D = new float[n2][n1];
				System.err.println("Reading image file.");
				input.read(img2D);

				/////////////////////////////////////////////////////////////////
				// Compute tensors
				System.err.println("Computing structure tensors.");
				et2D = ImageTensorFactory.getStructureTensors2D(lof,img2D);

				/////////////////////////////////////////////////////////////////
				// Write tensors file and close inpute
				System.err.println("Writing structure-tensors file.");
				output.setN(3,4);
				output.setOrigin(3,0.f);
				output.setDelta(3,0.f);
				output.write(et2D);
				output.close();
				input.close();
        break;
      case STRUCTURE_3D:
				System.err.println("RSFTensor: 3D STRUCTURE TENSORS");

				/////////////////////////////////////////////////////////////////
				// Read image file
				img3D = new float[n3][n2][n1];
				System.err.println("Reading image file.");
				input.read(img3D);

				/////////////////////////////////////////////////////////////////
				// Compute tensors
				System.err.println("Computing structure tensors.");
				et3D = ImageTensorFactory.getStructureTensors3D(lof,img3D);

				/////////////////////////////////////////////////////////////////
				// Write tensors file and close inpute
				System.err.println("Writing structure-tensors file.");
				output.setN(4,9);
				output.setOrigin(4,0.f);
				output.setDelta(4,0.f);
				output.write(et3D);
				output.close();
				input.close();
        break;
      case ORIENTED_2D:
				System.err.println("RSFTensor: 2D ORIENTED TENSORS");

				/////////////////////////////////////////////////////////////////
				// Read image file
				img2D = new float[n2][n1];
				System.err.println("Reading image file.");
				input.read(img2D);

				/////////////////////////////////////////////////////////////////
				// Compute tensors
				System.err.println("Computing oriented tensors.");
				et2D = ImageTensorFactory.getOrientedTensors2D(lof,img2D);

				/////////////////////////////////////////////////////////////////
				// Write tensors file and close inpute
				System.err.println("Writing oriented-tensors file.");
				output.setN(3,4);
				output.setOrigin(3,0.f);
				output.setDelta(3,0.f);
				output.write(et2D);
				output.close();
				input.close();
        break;
      case ORIENTED_3D:
				System.err.println("RSFTensor: 3D ORIENTED TENSORS");

				/////////////////////////////////////////////////////////////////
				// Read image file
				img3D = new float[n3][n2][n1];
				System.err.println("Reading image file.");
				input.read(img3D);

				/////////////////////////////////////////////////////////////////
				// Compute tensors
				System.err.println("Computing oriented tensors.");
				et3D = ImageTensorFactory.getOrientedTensors3D(lof,img3D);

				/////////////////////////////////////////////////////////////////
				// Write tensors file and close inpute
				System.err.println("Writing oriented-tensors file.");
				output.setN(4,9);
				output.setOrigin(4,0.f);
				output.setDelta(4,0.f);
				output.write(et3D);
				output.close();
				input.close();
        break;
      case EXCLUSION_2D:
				System.err.println("RSFTensor: 2D EXCLUSION TENSORS");
				/////////////////////////////////////////////////////////////////
				// get parameters
				float ecc2 = parms.getFloat("ecc",0.f);
				if(ecc2 <= 0.f) {
					System.err.println("An invalid eccentricity value was provided.");
					output.close();
					input.close();
					System.exit(-1);
				}

				/////////////////////////////////////////////////////////////////
				// Read tensors file
				img2D = new float[n2][n1];
				System.err.println("Reading image file.");
				input.read(img2D);

				/////////////////////////////////////////////////////////////////
				// Compute tensors
				System.err.println("Computing exclusion tensors.");
				et2D = ImageTensorFactory.getExclusionTensors2D(lof,img2D,ecc2);

				/////////////////////////////////////////////////////////////////
				// Write tensors file and close inpute
				System.err.println("Writing exclusion-tensors file.");
				output.setN(3,4);
				output.setOrigin(3,0.f);
				output.setDelta(3,0.f);
				output.write(et2D);
				output.close();
				input.close();
        break;
      case EXCLUSION_3D:
				System.err.println("RSFTensor: 3D EXCLUSION TENSORS");
				/////////////////////////////////////////////////////////////////
				// get parameters
				float ecc3 = parms.getFloat("ecc",0.f);
				if(ecc3 <= 0.f) {
					System.err.println("An invalid eccentricity value was provided.");
					output.close();
					input.close();
					System.exit(-1);
				}

				/////////////////////////////////////////////////////////////////
				// Read tensors file
				img3D = new float[n3][n2][n1];
				System.err.println("Reading image file.");
				input.read(img3D);

				/////////////////////////////////////////////////////////////////
				// Compute tensors
				System.err.println("Computing exclusion tensors.");
				et3D = ImageTensorFactory.getExclusionTensors3D(lof,img3D,ecc3);

				/////////////////////////////////////////////////////////////////
				// Write tensors file and close inpute
				System.err.println("Writing exclusion-tensors file.");
				output.setN(4,9);
				output.setOrigin(4,0.f);
				output.setDelta(4,0.f);
				output.write(et3D);
				output.close();
				input.close();
        break;
      case TENSOR_AT_COORD_2D:
				System.err.println("RSFTensor: 2D TENSORS-AT-COORDINATES");
				/////////////////////////////////////////////////////////////////
				// get parameters
				cfile = parms.getString("coordf","");
				if( cfile.equals("") ) {
					System.err.println("A coordinate file was not specified.");
					input.close();
					output.close();
					System.exit(-1);
				}

				/////////////////////////////////////////////////////////////////
				// Read tensors file
				et2D = new float[4][n2][n1];
				System.err.println("Reading tensors file.");
				input.read(et2D);


				/////////////////////////////////////////////////////////////////
				// Read coordinate file
				Input cInput2D = new Input("coordf");
				npicks = cInput2D.getN(1);
				if( npicks == 1 )
					System.err.println("Warning!! Only one coordinate was specifed.");
				ox1 = 0.f;
				ox2 = 0.f;
				dx1 = 0.f;
				dx2 = 0.f;
				ox1 = parms.getFloat("o1",0.f);
				ox2 = parms.getFloat("o2",0.f);
				ox3 = parms.getFloat("o3",0.f);
				dx1 = parms.getFloat("d1",0.f);
				dx2 = parms.getFloat("d2",0.f);
				coord = new float[2][npicks];
				System.err.println("Reading coordinate file: " + cfile);
				cInput2D.read(coord);
				float[][] etc = ImageTensorFactory.getTensorsByCoord2D(coord, et2D,
																															 ox1, ox2, 
																															 dx1, dx2);

				/////////////////////////////////////////////////////////////////
				// Write tensors at coordinates file and close inputs
				System.err.println("Writing tensors-at-coordinates file.");
				output.setN(1,npicks);
				output.setN(2,4);
				output.setN(3,1);
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);
				output.setDelta(1,0.f);
				output.setDelta(2,0.f);
				output.setDelta(3,1.f);
				output.write(etc);
				input.close();
				cInput2D.close();
				output.close();
        break;
      case TENSOR_AT_COORD_3D:
				System.err.println("RSFTensor: 3D TENSORS-AT-COORDINATES");
				/////////////////////////////////////////////////////////////////
				// get parameters
				cfile = parms.getString("coordf","");
				if( cfile.equals("") ) {
					System.err.println("A coordinate file was not specified.");
					input.close();
					output.close();
					System.exit(-1);
				}

				/////////////////////////////////////////////////////////////////
				// Read tensors file
				et3D = new float[9][n3][n2][n1];
				System.err.println("Reading tensors file.");
				input.read(et3D);


				/////////////////////////////////////////////////////////////////
				// Read coordinate file
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
				System.err.println("Reading coordinate file: " + cfile);
				cInput3D.read(coord);
				float[][] etc3D = ImageTensorFactory.getTensorsByCoord3D(coord, et3D,
																																 ox1, ox2, ox3, 
																																 dx1, dx2, dx3);

				/////////////////////////////////////////////////////////////////
				// Write tensors at coordinates file and close inputs
				System.err.println("Writing tensors-at-coordinates file.");
				output.setN(1,npicks);
				output.setN(2,9);
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
				output.write(etc3D);
				input.close();
				cInput3D.close();
				output.close();
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
			if(sopt_str.equals("metric"))
				return SwitchOpt.METRIC_3D;
			if(sopt_str.equals("struc"))
				return SwitchOpt.STRUCTURE_3D;
			if(sopt_str.equals("orien"))
				return SwitchOpt.ORIENTED_3D;
			if(sopt_str.equals("excld"))
				return SwitchOpt.EXCLUSION_3D;
			if(sopt_str.equals("coord"))
				return SwitchOpt.TENSOR_AT_COORD_3D;
			return SwitchOpt.DEFAULT;
		} 
		if(sopt_str.equals("metric"))
			return SwitchOpt.METRIC_2D;
		if(sopt_str.equals("struc"))
			return SwitchOpt.STRUCTURE_2D;
		if(sopt_str.equals("orien"))
			return SwitchOpt.ORIENTED_2D;
		if(sopt_str.equals("excld"))
			return SwitchOpt.EXCLUSION_2D;
		if(sopt_str.equals("coord"))
			return SwitchOpt.TENSOR_AT_COORD_2D;
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

}


