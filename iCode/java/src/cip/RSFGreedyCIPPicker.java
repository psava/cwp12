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

import edu.mines.jtk.dsp.*;

import rsf.RSF;
import rsf.Input;
import rsf.Output;

/**
 * Picks CIP locations using a greedy picking hueristic. Pixels or voxels 
 * are sorted by nondecreasing order of priority. A mapping of each
 * pixel's/voxel's index is maintained. The pixel/voxel with the highest
 * priority is picked first and the pixels/voxels around the pick (within
 * an exclusion zone) are then excluded from being picked.  There are two
 * Types of exclusion zones: grid-aligned ellipses/ellipsoids or structure
 * aligned noneuclidean ellipses/ellipsoids as defined by an eikonal
 * equation used in TimeSolver2/TimeSolver3 from the Mines JTK.
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.11.14
 */
public class RSFGreedyCIPPicker
{

	static { 
	  System.loadLibrary("jrsf");
	}

  enum SwitchOpt{
		INDEX_2D,
		INDEX_3D,
		COORD_2D,
		COORD_3D,
		TENSOR_COORD_2D,
		TENSOR_COORD_3D,
    DEFAULT
  };

	public static void main (String[] args){
		RSF parms = new RSF(args);
		Input input = new Input("in");
		Output output = new Output("out");

    String opt_str = parms.getString("opt","");
    if(opt_str.equals("")) {
			System.err.println("A picker option was not provided.");
			output.close();
			input.close();
			System.exit(-1);
		}

    String etype_str = parms.getString("etype","none");

		float r1,r2,r3;
    r1 = parms.getFloat("r1",0.f);
    if(r1 <= 0.f) {
			System.err.println("r1 must be greater than zero.");
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


    SwitchOpt switch_opt = getSwitchOpt(opt_str,etype_str,dim);


		//////////////////////////////////////////////////////
		// Get dims
		int n1,n2,n3;
		n1 = input.getN(1);
		n2 = input.getN(2);
		n3 = input.getN(3);
		

		//////////////////////////////////////////////////////
		// Get deltas
		float dx1,dx2,dx3;
		dx1 = input.getDelta(1);
		dx2 = input.getDelta(2);
		dx3 = input.getDelta(3);


		//////////////////////////////////////////////////////
		// Get origins
		float ox1,ox2,ox3;
    ox1 = input.getOrigin(1);
    ox2 = input.getOrigin(2);
		ox3 = input.getOrigin(3);


		//////////////////////////////////////////////////////////
		// Carryout desired option
		int npicks = 0;
		float[][]   prob2D;
		float[][][] prob3D;
		float[][]   cipCoord;
		int[][]   cipIndex;
		String tfile;
		GreedyCIPPicker2 gCIPPicker2;
		GreedyCIPPicker3 gCIPPicker3;
    switch (switch_opt)
    {
      case INDEX_2D:
				//FIXME: SWIG RSF interface does not handel integer data files.
				/*
				prob2D = new float[n2][n1];
				input.read(prob2D);
				gCIPPicker2 = new GreedyCIPPicker2(prob2D, radius);
				cipIndex = gCIPPicker2.getCIPIndices();
				npicks = gCIPPicker2.getNumPicks();
				output.setN(1,npicks);
				output.setN(2,2);
				output.setN(3,1);//rsf problem with 2D images keeping 3D data info
				output.setDelta(1,0.f);
				output.setDelta(2,0.f);
				output.setDelta(3,0.f);//rsf problem with 2D images keeping 3D data info
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);//rsf problem with 2D images keeping 3D data info
				output.write(cipIndex);
				output.close();
				input.close();
				*/
        break;
      case INDEX_3D:
				//FIXME: SWIG RSF interface does not handel integer data files.
				/*
				prob3D = new float[n3][n2][n1];
				input.read(prob3D);
				gCIPPicker3 = new GreedyCIPPicker3(prob3D, radius);
				cipIndex = gCIPPicker3.getCIPIndices();
				npicks = gCIPPicker3.getNumPicks();
				output.setN(1,npicks);
				output.setN(2,3);
				output.setN(3,1);
				output.setDelta(1,0.f);
				output.setDelta(2,0.f);
				output.setDelta(3,0.f);
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);
				output.write(cipIndex);
				output.close();
				input.close();
				*/
        break;
      case COORD_2D:
				System.err.println("RSFGreedyCIPPicker: 2D COORDINATES");
				/////////////////////////////////////////////////////////
				// Get parameters
				r2 = parms.getFloat("r2",0.f);
				if(r2 <= 0.f) {
					System.err.println("r2 must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}
				if( dx1 <= 0.f || dx2 <= 0.f ) {
					System.err.println("The dx1, dx2 options must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}

				/////////////////////////////////////////////////////////
				// Read priority map
				prob2D = new float[n2][n1];
				System.err.println("Reading priority map.");
				input.read(prob2D);

				/////////////////////////////////////////////////////////
				// Get Picks
				System.err.println("Picking CIP locations.");
				System.err.println("r1 = " + (int)(r1+0.5f));
				System.err.println("r2 = " + (int)(r2+0.5f));
				EllipseCIPExcluder eExcluder2 = new EllipseCIPExcluder((int)(r1+0.5f),
																															 (int)(r2+0.5f),
																															 n1,n2,
																															 false);
				gCIPPicker2 = new GreedyCIPPicker2(prob2D,eExcluder2);

				/////////////////////////////////////////////////////////
				// Computing Coordinates
				System.err.println("Picking CIP locations.");
				cipCoord = gCIPPicker2.getCIPCoord(ox1,ox2,dx1,dx2);
				npicks = gCIPPicker2.getNumPicks();
				System.err.println("Number of picks = " + npicks + ".");

				/////////////////////////////////////////////////////////
				// Write Coordinates and close inputs
				System.err.println("Writing 2D Coordinates");
				output.setN(1,npicks);
				output.setN(2,2);
				output.setN(3,1);//rsf problem with 2D images keeping 3D data info
				output.setDelta(1,0.f);
				output.setDelta(2,0.f);
				output.setDelta(3,0.f);//rsf problem with 2D images keeping 3D data info
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);//rsf problem with 2D images keeping 3D data info
				output.write(cipCoord);
				output.close();
				input.close();
        break;
      case COORD_3D:
				System.err.println("RSFGreedyCIPPicker: 3D COORDINATES");
				/////////////////////////////////////////////////////////
				// Get parameters
				r2 = parms.getFloat("r2",0.f);
				if(r2 <= 0.f) {
					System.err.println("r2 must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}
				r3 = parms.getFloat("r3",0.f);
				if(r3 <= 0.f) {
					System.err.println("r3 must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}
				if( dx1 <= 0.f || dx2 <= 0.f || dx3 <= 0.f ) {
					System.err.println("The dx1, dx2, dx3 options must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}

				/////////////////////////////////////////////////////////
				// Read priority map
				prob3D = new float[n3][n2][n1];
				System.err.println("Reading priority map.");
				input.read(prob3D);

				/////////////////////////////////////////////////////////
				// Get Picks
				System.err.println("Picking CIP locations.");
				EllipsoidCIPExcluder eExcluder3 = new EllipsoidCIPExcluder((int)(r1+0.5f),
																																	 (int)(r2+0.5f),
																																	 (int)(r3+0.5f),
																																	 n1,n2,n3,
																																	 false);
				gCIPPicker3 = new GreedyCIPPicker3(prob3D,eExcluder3);

				/////////////////////////////////////////////////////////
				// Computing Coordinates
				System.err.println("Picking CIP locations.");
				cipCoord = gCIPPicker3.getCIPCoord(ox1,ox2,ox3,dx1,dx2,dx3);
				npicks = gCIPPicker3.getNumPicks();
				System.err.println("Number of picks = " + npicks + ".");

				/////////////////////////////////////////////////////////
				// Write Coordinates and close inputs
				output.setN(1,npicks);
				output.setN(2,3);
				output.setN(3,1);
				output.setDelta(1,0.f);
				output.setDelta(2,0.f);
				output.setDelta(3,0.f);
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);
				output.write(cipCoord);
				output.close();
				input.close();
        break;
      case TENSOR_COORD_2D:
				System.err.println("RSFGreedyCIPPicker: 2D TENSOR-COORDINATES");
				/////////////////////////////////////////////////////////
				// Get parameters
				tfile = parms.getString("tfile","");
				if(tfile.equals("")) {
					System.err.println("A tensor file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				if( dx1 <= 0.f || dx2 <= 0.f) {
					System.err.println("The dx1 and dx2 options must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}

				/////////////////////////////////////////////////////////
				// Read priority map
				prob2D = new float[n2][n1];
				float[][][] et = new float[4][n2][n1];
				System.err.println("Reading priority map.");
				input.read(prob2D);

				/////////////////////////////////////////////////////////
				// Read tensor file
				Input tenInput = new Input("tfile");
				System.err.println("Reading tensor file: " + tfile);
				tenInput.read(et);

				/////////////////////////////////////////////////////////
				// Load tensors
				System.err.println("Loading tensors.");
				EigenTensors2 et2 = new EigenTensors2(n1,n2);
				float[] u = new float[2];
				for(int i2=0;i2<n2;++i2)
					for(int i1=0;i1<n1;++i1) {
						u[0] = et[0][i2][i1];
						u[1] = et[1][i2][i1];
						et2.setEigenvectorU(i1,i2,u);
					}
				et2.setEigenvalues(et[2],et[3]);


				/////////////////////////////////////////////////////////
				// Write Coordinates and close inputs
				System.err.println("Writing 2D Coordinates");
				TensorEllipseCIPExcluder teExcluder2 = new TensorEllipseCIPExcluder(et2,r1,n1,n2,false);
				gCIPPicker2 = new GreedyCIPPicker2(prob2D,teExcluder2);
				cipCoord = gCIPPicker2.getCIPCoord(ox1,ox2,dx1,dx2);
				npicks = gCIPPicker2.getNumPicks();
				output.setN(1,npicks);
				output.setN(2,2);
				output.setN(3,1);//rsf problem with 2D images keeping 3D data info
				output.setDelta(1,0.f);
				output.setDelta(2,0.f);
				output.setDelta(3,0.f);//rsf problem with 2D images keeping 3D data info
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);//rsf problem with 2D images keeping 3D data info
				output.write(cipCoord);
				output.close();
				input.close();
        break;
      case TENSOR_COORD_3D:
				System.err.println("RSFGreedyCIPPicker: 3D TENSOR-COORDINATES");
				/////////////////////////////////////////////////////////
				// Get parameters
				tfile = parms.getString("tfile","");
				if(tfile.equals("")) {
					System.err.println("A tensor file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				if( dx1 <= 0.f || dx2 <= 0.f || dx3 <= 0.f ) {
					System.err.println("The dx1, dx2, dx3 options must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}

				////////////////////////////////////////
				// read priority map
				prob3D = new float[n3][n2][n1];
				float[][][][] rsfet = new float[9][n3][n2][n1];
				System.err.println("Reading priority file.");
				input.read(prob3D);

				////////////////////////////////////////
				// read tensors
				System.err.println("Reading tensor file.");
				Input tenInput3D = new Input("tfile");
				tenInput3D.read(rsfet);

				////////////////////////////////////////
				// load tensors
				System.err.println("Loading tensors.");
				EigenTensors3 et3 = new EigenTensors3(n1,n2,n3,true);
				float u1,u2,u3,w1,w2,w3;
				for(int i3=0;i3<n3;++i3)
					for(int i2=0;i2<n2;++i2)
						for(int i1=0;i1<n1;++i1) {
						  u1 = rsfet[0][i3][i2][i1];
						  u2 = rsfet[1][i3][i2][i1];
						  u3 = rsfet[2][i3][i2][i1];
						  w1 = rsfet[3][i3][i2][i1];
						  w2 = rsfet[4][i3][i2][i1];
						  w3 = rsfet[5][i3][i2][i1];
						  et3.setEigenvectorU(i1,i2,i3,u1,u2,u3);
						  et3.setEigenvectorW(i1,i2,i3,w1,w2,w3);
						}
				et3.setEigenvalues(rsfet[6],rsfet[7],rsfet[8]);

				////////////////////////////////////////
				// Get picks
				System.err.println("Picking CIP locations.");
				TensorEllipsoidCIPExcluder teExcluder3 = new TensorEllipsoidCIPExcluder(et3,r1,n1,n2,n3,false);
				gCIPPicker3 = new GreedyCIPPicker3(prob3D,teExcluder3);
				cipCoord = gCIPPicker3.getCIPCoord(ox1,ox2,ox3,dx1,dx2,dx3);
				npicks = gCIPPicker3.getNumPicks();
				System.err.println("Number of picks = " + npicks);

				////////////////////////////////////////
				// write and close inputs 
				System.err.println("Writing pick-coordinates file.");
				output.setN(1,npicks);
				output.setN(2,3);
				output.setN(3,1);
				output.setDelta(1,0.f);
				output.setDelta(2,0.f);
				output.setDelta(3,0.f);
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);
				output.write(cipCoord);
				output.close();
				input.close();
				tenInput3D.close();
        break;
      default:
				System.err.println("Cannot find a match for the picker option provided.");
				output.close();
				input.close();
				System.exit(-1);
    }

  } //end main()


	/**********************************************************************/
	// 
	// Private Member Functions
	// 
	/**********************************************************************/

	//////////////////////////////////////////////////////
	// Translate argument option into an ENUM
  private static SwitchOpt getSwitchOpt(String sopt_str, String etype_str, int dim) {
		if(dim == 3) {
			if(sopt_str.equals("coord"))
				if(etype_str.equals("tensor"))
					return SwitchOpt.TENSOR_COORD_3D;
				else
					return SwitchOpt.COORD_3D;
			if(sopt_str.equals("index"))
				return SwitchOpt.INDEX_3D;
			return SwitchOpt.DEFAULT;
		}
		if(sopt_str.equals("coord"))
			if(etype_str.equals("tensor"))
				return SwitchOpt.TENSOR_COORD_2D;
			else
				return SwitchOpt.COORD_2D;
		if(sopt_str.equals("index"))
			return SwitchOpt.INDEX_2D;
		return SwitchOpt.DEFAULT;
	} //end getSwitchOpt()

} //end class
