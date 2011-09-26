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

import cip.util.*;

/**
 * Computes the exclusion-zone semblance for CIP picks.
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.11.14
 */
public class RSFEZSemblanceSort
{

	static { 
	  System.loadLibrary("jrsf");
	}

  enum SwitchOpt{
		ELLP_EZSEM_2D,
		ELLP_EZSEM_3D,
		TENSOR_EZSEM_2D,
		TENSOR_EZSEM_3D,
    DEFAULT
  };

	public static void main (String[] args){
		RSF parms = new RSF(args);
		Input input = new Input("in");
		Output output = new Output("out");

    String opt_str = parms.getString("opt","");
    if(opt_str.equals("")) {
			System.err.println("An exclusion-zone-semblance option was not provided.");
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
		PointOfInterest[] poi;
		PointOfInterest[] ezpoi;
		int[] indices;
		float[][]   coord;
		float[][]   pri2D;
		float[][][] pri3D;
		float rmax = 0.f;
		int r1 = 0;
		int r2 = 0;
		int r3 = 0;
		int f = 0;
		int j = 1;
		int npicks = 0;
		String tfile,cfile;
    switch (switch_opt)
    {
      case ELLP_EZSEM_2D:
				System.err.println("RSFEZSemblanceSort: ELLIPSE SEMBLANCE");
				/////////////////////////////////////////////////
				// get parameters
				cfile = parms.getString("cfile","");
				if(cfile.equals("")) {
					System.err.println("A coordinate file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				r1 = (int)(parms.getFloat("r1",0.f) + 0.5f);
				if(r1 <= 0) {
					System.err.println("r1 must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}
				r2 = (int)(parms.getFloat("r2",0.f) + 0.5f);
				if(r2 <= 0) {
					System.err.println("r2 must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}

				/////////////////////////////////////////////////
				// read priority file
				pri2D = new float[n2][n1];
				System.err.println("Reading priority-map file");
				input.read(pri2D);

				/////////////////////////////////////////////////
				// read coordinate file
				Input cInput = new Input("cfile");
				npicks = cInput.getN(1);
				if(npicks <= 0) {
					System.err.println("The coordinate file indicates that there were zero picks");
					output.close();
					input.close();
					System.exit(-1);
				}
				coord = new float[2][npicks];
				System.err.println("Reading coordinate file: " + cfile);
				cInput.read(coord);

				/////////////////////////////////////////////////
				// create poi array
				System.err.println("Loading POI array.");
				poi = new PointOfInterest[npicks];
				indices = new int[2];
				for(int i=0; i<npicks; ++i) {
						indices[0] = (int)((coord[0][i]-ox1)/dx1 + 0.5f);
						indices[1] = (int)((coord[1][i]-ox2)/dx2 + 0.5f);
						float val = pri2D[indices[1]][indices[0]];
						PointOfInterest tpoi = new PointOfInterest(val,2,indices);
						poi[i] = tpoi;
				}

				/////////////////////////////////////////////////////////
				// calculate exclusion-zone-semblance
				System.err.println("Calculating exclusion-zone semblances.");
				EllipseCIPEZSemblance ezSemblance = new EllipseCIPEZSemblance(r1,r2,n1,n2);
				EZSemblanceSort.sortByEZSemblance(pri2D,poi,ezSemblance);

				/////////////////////////////////////////////////////////
				// load resorted coordinates
				System.err.println("Calculating resorted coordinates.");
				loadResortedCoords(coord,poi,dx1,dx2,ox1,ox2);

				/////////////////////////////////////////////////////////
				// Write image file and close inputs
				System.err.println("Writing resorted-pick-coordinates file.");
				output.setN(1,npicks);
				output.setN(2,2);
				output.setN(3,1);//rsf problem with 2D data keeping 3D data info
				output.setDelta(1,1.f);
				output.setDelta(2,1.f);
				output.setDelta(3,1.f);//rsf problem with 2D data keeping 3D data info
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);//rsf problem with 2D data keeping 3D data info
				output.write(coord);
				output.close();
				input.close();
				cInput.close();
        break;
      case ELLP_EZSEM_3D:
				System.err.println("RSFEZSemblanceSort: ELLIPSOID SEMBLANCE");
				/////////////////////////////////////////////////
				// get parameters
				cfile = parms.getString("cfile","");
				if(cfile.equals("")) {
					System.err.println("A coordinate file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				r1 = (int)(parms.getFloat("r1",0.f) + 0.5f);
				if(r1 <= 0) {
					System.err.println("r1 must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}
				r2 = (int)(parms.getFloat("r2",0.f) + 0.5f);
				if(r2 <= 0) {
					System.err.println("r2 must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}
				r3 = (int)(parms.getFloat("r3",0.f) + 0.5f);
				if(r3 <= 0) {
					System.err.println("r3 must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}

				/////////////////////////////////////////////////
				// read priority file
				pri3D = new float[n3][n2][n1];
				System.err.println("Reading priority-map file");
				input.read(pri3D);

				/////////////////////////////////////////////////
				// read coordinate file
				Input cInput3D = new Input("cfile");
				npicks = cInput3D.getN(1);
				if(npicks <= 0) {
					System.err.println("The coordinate file indicates that there were zero picks");
					output.close();
					input.close();
					System.exit(-1);
				}
				coord = new float[3][npicks];
				System.err.println("Reading coordinate file: " + cfile);
				cInput3D.read(coord);

				/////////////////////////////////////////////////
				// create poi array
				System.err.println("Loading POI array.");
				poi = new PointOfInterest[npicks];
				indices = new int[3];
				for(int i=0; i<npicks; ++i) {
						indices[0] = (int)((coord[0][i]-ox1)/dx1 + 0.5f);
						indices[1] = (int)((coord[1][i]-ox2)/dx2 + 0.5f);
						indices[2] = (int)((coord[2][i]-ox3)/dx3 + 0.5f);
						float val = pri3D[indices[2]][indices[1]][indices[0]];
						PointOfInterest tpoi = new PointOfInterest(val,3,indices);
						poi[i] = tpoi;
				}

				/////////////////////////////////////////////////////////
				// calculate exclusion-zone-semblance
				System.err.println("Calculating exclusion-zone semblances.");
				EllipsoidCIPEZSemblance ezSemblance3 = new EllipsoidCIPEZSemblance(r1,r2,r3,
																																						n1,n2,n3);
				EZSemblanceSort.sortByEZSemblance(pri3D,poi,ezSemblance3);

				/////////////////////////////////////////////////////////
				// load resorted coordinates
				System.err.println("Calculating resorted coordinates.");
				loadResortedCoords(coord,poi,dx1,dx2,dx3,ox1,ox2,ox3);

				/////////////////////////////////////////////////////////
				// Write image file and close inputs
				System.err.println("Writing resorted-pick-coordinates file.");
				output.setN(1,npicks);
				output.setN(2,3);
				output.setN(3,1);
				output.setDelta(1,1.f);
				output.setDelta(2,1.f);
				output.setDelta(3,1.f);
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);
				output.write(coord);
				output.close();
				input.close();
				cInput3D.close();
        break;
      case TENSOR_EZSEM_2D:
				System.err.println("RSFEZSemblanceSort: TENSOR-ELLIPSE SEMBLANCE");
				/////////////////////////////////////////////////
				// get parameters
				cfile = parms.getString("cfile","");
				if(cfile.equals("")) {
					System.err.println("A coordinate file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				tfile = parms.getString("tfile","");
				if(tfile.equals("")) {
					System.err.println("A tensor file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				rmax = parms.getFloat("rmax",0.f);
				if(rmax <= 0.f) {
					System.err.println("rmax must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}

				/////////////////////////////////////////////////
				// read coordinate file
				Input cInput2 = new Input("cfile");
				npicks = cInput2.getN(1);
				if(npicks <= 0) {
					System.err.println("The coordinate file indicates that there were zero picks");
					output.close();
					input.close();
					System.exit(-1);
				}
				coord = new float[2][npicks];
				System.err.println("Reading coordinate file: " + cfile);
				cInput2.read(coord);

				/////////////////////////////////////////////////
				// read priority file
				pri2D = new float[n2][n1];
				System.err.println("Reading priority-map file");
				input.read(pri2D);

				/////////////////////////////////////////////////
				// read tensors file
				float[][][] et = new float[4][n2][n1];
				Input tenInput = new Input("tfile");
				System.err.println("Reading tensor file: " + tfile);
				tenInput.read(et);

				/////////////////////////////////////////////////
				// create poi array
				System.err.println("Loading POI array.");
				poi = new PointOfInterest[npicks];
				indices = new int[2];
				for(int i=0; i<npicks; ++i) {
						indices[0] = (int)((coord[0][i]-ox1)/dx1 + 0.5f);
						indices[1] = (int)((coord[1][i]-ox2)/dx2 + 0.5f);
						float val = pri2D[indices[1]][indices[0]];
						PointOfInterest tpoi = new PointOfInterest(val,2,indices);
						poi[i] = tpoi;
				}

				/////////////////////////////////////////////////
				// create tensors
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
				// calculate exclusion-zone-semblance
				System.err.println("Calculating exclusion-zone semblances.");
				TensorEllipseCIPEZSemblance tezSemblance = new TensorEllipseCIPEZSemblance(et2,rmax,n1,n2);
				EZSemblanceSort.sortByEZSemblance(pri2D,poi,tezSemblance);

				/////////////////////////////////////////////////////////
				// load resorted coordinates
				System.err.println("Calculating resorted coordinates.");
				loadResortedCoords(coord,poi,dx1,dx2,ox1,ox2);

				/////////////////////////////////////////////////////////
				// Write coord file and close inputs
				System.err.println("Writing resorted-pick-coordinates file.");
				output.setN(1,npicks);
				output.setN(2,2);
				output.setN(3,1);//rsf problem with 2D data keeping 3D data info
				output.setDelta(1,1.f);
				output.setDelta(2,1.f);
				output.setDelta(3,1.f);//rsf problem with 2D data keeping 3D data info
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);//rsf problem with 2D data keeping 3D data info
				output.write(coord);
				output.close();
				input.close();
				cInput2.close();
				tenInput.close();
        break;
      case TENSOR_EZSEM_3D:
				System.err.println("RSFEZSemblanceSort: TENSOR-ELLIPSOID SEMBLANCE");
				/////////////////////////////////////////////////
				// get parameters
				cfile = parms.getString("cfile","");
				if(cfile.equals("")) {
					System.err.println("A coordinate file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				tfile = parms.getString("tfile","");
				if(tfile.equals("")) {
					System.err.println("A tensor file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				rmax = parms.getFloat("rmax",0.f);
				if(rmax <= 0.f) {
					System.err.println("rmax must be greater than zero.");
					output.close();
					input.close();
					System.exit(-1);
				}

				/////////////////////////////////////////////////
				// read priority file
				pri3D = new float[n3][n2][n1];
				System.err.println("Reading priority-map file");
				input.read(pri3D);

				/////////////////////////////////////////////////
				// read coordinate file
				Input tcInput3D = new Input("cfile");
				npicks = tcInput3D.getN(1);
				if(npicks <= 0) {
					System.err.println("The coordinate file indicates that there were zero picks");
					output.close();
					input.close();
					System.exit(-1);
				}
				coord = new float[3][npicks];
				System.err.println("Reading coordinate file: " + cfile);
				tcInput3D.read(coord);

				/////////////////////////////////////////////////
				// read tensors file
				float[][][][] et3D = new float[9][n3][n2][n1];
				System.err.println("Reading tensor file: " + tfile);
				Input tenInput3D = new Input("tfile");
				tenInput3D.read(et3D);

				/////////////////////////////////////////////////
				// create poi array
				System.err.println("Loading POI array.");
				poi = new PointOfInterest[npicks];
				indices = new int[3];
				for(int i=0; i<npicks; ++i) {
						indices[0] = (int)((coord[0][i]-ox1)/dx1 + 0.5f);
						indices[1] = (int)((coord[1][i]-ox2)/dx2 + 0.5f);
						indices[2] = (int)((coord[2][i]-ox3)/dx3 + 0.5f);
						float val = pri3D[indices[2]][indices[1]][indices[0]];
						PointOfInterest tpoi = new PointOfInterest(val,3,indices);
						poi[i] = tpoi;
				}

				/////////////////////////////////////////////////
				// create tensors
				System.err.println("Loading tensors.");
				EigenTensors3 et3 = new EigenTensors3(n1,n2,n3,true);
				float u1,u2,u3,w1,w2,w3;
				for(int i3=0;i3<n3;++i3)
					for(int i2=0;i2<n2;++i2)
						for(int i1=0;i1<n1;++i1) {
						  u1 = et3D[0][i3][i2][i1];
						  u2 = et3D[1][i3][i2][i1];
						  u3 = et3D[2][i3][i2][i1];
						  w1 = et3D[3][i3][i2][i1];
						  w2 = et3D[4][i3][i2][i1];
						  w3 = et3D[5][i3][i2][i1];
						  et3.setEigenvectorU(i1,i2,i3,u1,u2,u3);
						  et3.setEigenvectorW(i1,i2,i3,w1,w2,w3);
						}
				et3.setEigenvalues(et3D[6],et3D[7],et3D[8]);

				/////////////////////////////////////////////////////////
				// calculate exclusion-zone-semblance
				System.err.println("Calculating exclusion-zone semblances.");
				TensorEllipsoidCIPEZSemblance tezSemblance3 = 
																			new TensorEllipsoidCIPEZSemblance(et3,rmax,n1,n2,n3);
				EZSemblanceSort.sortByEZSemblance(pri3D,poi,tezSemblance3);

				/////////////////////////////////////////////////////////
				// load resorted coordinates
				System.err.println("Calculating resorted coordinates.");
				loadResortedCoords(coord,poi,dx1,dx2,dx3,ox1,ox2,ox3);

				/////////////////////////////////////////////////////////
				// Write coordinate file and close inputs
				System.err.println("Writing image with exclusion zones.");
				output.setN(1,npicks);
				output.setN(2,3);
				output.setN(3,1);
				output.setDelta(1,1.f);
				output.setDelta(2,1.f);
				output.setDelta(3,1.f);
				output.setOrigin(1,0.f);
				output.setOrigin(2,0.f);
				output.setOrigin(3,0.f);
				output.write(coord);
				output.close();
				input.close();
				tcInput3D.close();
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
  private static SwitchOpt getSwitchOpt(String sopt_str, int dim) {
		if(dim == 3) {
			if(sopt_str.equals("tensor"))
				return SwitchOpt.TENSOR_EZSEM_3D;
			if(sopt_str.equals("ellp"))
				return SwitchOpt.ELLP_EZSEM_3D;
			return SwitchOpt.DEFAULT;
		}
		if(sopt_str.equals("tensor"))
			return SwitchOpt.TENSOR_EZSEM_2D;
		if(sopt_str.equals("ellp"))
			return SwitchOpt.ELLP_EZSEM_2D;
		return SwitchOpt.DEFAULT;
	} //end getSwitchOpt()


	/////////////////////////////////////////////////////////
	// load resorted coordinates
	private static void loadResortedCoords(float[][] coord, 
																				 PointOfInterest[] ezpoi,
																				 float dx1, float dx2, 
																				 float ox1, float ox2) {
		for(int i=0; i<ezpoi.length; ++i) {
			coord[0][i] = ox1 + dx1*ezpoi[i].getIndex(0);
			coord[1][i] = ox2 + dx2*ezpoi[i].getIndex(1);
		}
	}

	/////////////////////////////////////////////////////////
	// load resorted coordinates
	private static void loadResortedCoords(float[][] coord, 
																				 PointOfInterest[] ezpoi,
																				 float dx1, float dx2, float dx3, 
																				 float ox1, float ox2, float ox3) {
		System.err.println("coord[0].length = " + coord[0].length);
		System.err.println("ezpoi.length = " + ezpoi.length);
		System.err.println("ox1,ox2,ox3 = " +ox1+","+ox2+","+ox3);
		System.err.println("dx1,dx2,dx3 = " +dx1+","+dx2+","+dx3);
		for(int i=0; i<ezpoi.length; ++i) {
			//System.err.println("i = " + i);
			coord[0][i] = ox1 + dx1*ezpoi[i].getIndex(0);
			coord[1][i] = ox2 + dx2*ezpoi[i].getIndex(1);
			coord[2][i] = ox3 + dx3*ezpoi[i].getIndex(2);
		}
	}

} //end class
