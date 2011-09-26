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
 * Marks Exclusion zones at pick locations.
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.11.14
 */
public class RSFExclusionZones
{

	static { 
	  System.loadLibrary("jrsf");
	}

  enum SwitchOpt{
		ELLP_EXCLD_2D_ONE,
		ELLP_EXCLD_3D_ONE,
		TENSOR_EXCLD_2D_ONE,
		TENSOR_EXCLD_3D_ONE,
		ELLP_EXCLD_2D,
		ELLP_EXCLD_3D,
		TENSOR_EXCLD_2D,
		TENSOR_EXCLD_3D,
    DEFAULT
  };

	public static void main (String[] args){
		RSF parms = new RSF(args);
		Input input = new Input("in");
		Output output = new Output("out");

    String opt_str = parms.getString("opt","");
    if(opt_str.equals("")) {
			System.err.println("An exclusion-zone option was not provided.");
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
		PointOfInterest[] wpoi;
		int[] indices;
		float[][]   coord;
		float[][]   pri2D;
		float[][][] pri3D;
		float[][]   image2D;
		float[][][] image3D;
		float[][]   zimage2D;
		float[][][] zimage3D;
		float rmax = 0.f;
		int r1 = 0;
		int r2 = 0;
		int r3 = 0;
		int f = 0;
		int j = 1;
		int npicks = 0;
		float c1 = 0.f;
		float c2 = 0.f;
		float c3 = 0.f;
		float sval = 0.f;
		String tfile,pfile,cfile;
		GreedyCIPPicker2 gCIPPicker2;
		GreedyCIPPicker3 gCIPPicker3;
    switch (switch_opt)
    {
      case ELLP_EXCLD_2D_ONE:
				System.err.println("RSFExclusionZones: SIGNLE ELLIPSE EXCLUDER");
				/////////////////////////////////////////////////
				// get parameters
				pfile = parms.getString("pfile","");
				if(pfile.equals("")) {
					System.err.println("A priority-map file was not specified.");
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
				c1 = parms.getFloat("c1",0.f);
				System.err.println("c1 = " + c1);
				c2 = parms.getFloat("c2",0.f);
				System.err.println("c2 = " + c2);


				/////////////////////////////////////////////////
				// read image file to mark with zones
				image2D = new float[n2][n1];
				System.err.println("Reading image to be marked.");
				input.read(image2D);


				/////////////////////////////////////////////////
				// read priority file
				pri2D = new float[n2][n1];
				Input spriInput = new Input("pfile");
				System.err.println("Reading priority-map file: " + pfile);
				spriInput.read(pri2D);


				/////////////////////////////////////////////////
				// create poi array
				poi = new PointOfInterest[1];
				indices = new int[2];
				indices[0] = (int)((c1-ox1)/dx1 + 0.5f);
				indices[1] = (int)((c2-ox2)/dx2 + 0.5f);
				sval = pri2D[indices[1]][indices[0]];
				poi[0] = new PointOfInterest(sval,2,indices);


				/////////////////////////////////////////////////////////
				// window picks and create image with exclusion zone
				System.err.println("Marking exclusion zone.");
				EllipseCIPExcluder seExcluder2 = new EllipseCIPExcluder(r1,r2,n1,n2,false);
				zimage2D = ExclusionZones.showZones(image2D,poi,seExcluder2);


				/////////////////////////////////////////////////////////
				// Write image file and close inputs
				System.err.println("Writing image marked with exclusion zone.");
				output.setN(1,n1);
				output.setN(2,n2);
				output.setN(3,1);//rsf problem with 2D images keeping 3D data info
				output.setDelta(1,dx1);
				output.setDelta(2,dx2);
				output.setDelta(3,1.f);//rsf problem with 2D images keeping 3D data info
				output.setOrigin(1,ox1);
				output.setOrigin(2,ox2);
				output.setOrigin(3,0.f);//rsf problem with 2D images keeping 3D data info
				output.write(zimage2D);
				output.close();
				input.close();
				spriInput.close();
        break;
      case ELLP_EXCLD_3D_ONE:
				System.err.println("RSFExclusionZones: SINGLE ELLIPSOID EXCLUDER");
				/////////////////////////////////////////////////
				// get parameters
				pfile = parms.getString("pfile","");
				if(pfile.equals("")) {
					System.err.println("A priority-map file was not specified.");
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
				c1 = parms.getFloat("c1",0.f);
				System.err.println("c1 = " + c1);
				c2 = parms.getFloat("c2",0.f);
				System.err.println("c2 = " + c2);
				c3 = parms.getFloat("c3",0.f);
				System.err.println("c3 = " + c3);


				/////////////////////////////////////////////////
				// read image file to mark with zones
				image3D = new float[n3][n2][n1];
				System.err.println("Reading image to be marked.");
				input.read(image3D);


				/////////////////////////////////////////////////
				// read priority file
				pri3D = new float[n3][n2][n1];
				Input spriInput3D = new Input("pfile");
				System.err.println("Reading priority-map file: " + pfile);
				spriInput3D.read(pri3D);


				/////////////////////////////////////////////////
				// create poi array
				poi = new PointOfInterest[1];
				indices = new int[3];
				indices[0] = (int)((c1-ox1)/dx1 + 0.5f);
				indices[1] = (int)((c2-ox2)/dx2 + 0.5f);
				indices[2] = (int)((c3-ox3)/dx3 + 0.5f);
				sval = pri3D[indices[2]][indices[1]][indices[0]];
				poi[0] = new PointOfInterest(sval,3,indices);


				/////////////////////////////////////////////////////////
				// window picks and create image with exclusion zones
				System.err.println("Marking exclusion zone.");
				EllipsoidCIPExcluder seExcluder3 = new EllipsoidCIPExcluder(r1,r2,r3,n1,n2,n3,false);
				zimage3D = ExclusionZones.showZones(image3D,poi,seExcluder3);


				/////////////////////////////////////////////////////////
				// Write image file and close inputs
				System.err.println("Writing image marked with exclusion zone.");
				output.setN(1,n1);
				output.setN(2,n2);
				output.setN(3,n3);
				output.setDelta(1,dx1);
				output.setDelta(2,dx2);
				output.setDelta(3,dx3);
				output.setOrigin(1,ox1);
				output.setOrigin(2,ox2);
				output.setOrigin(3,ox3);
				output.write(zimage3D);
				output.close();
				input.close();
				spriInput3D.close();
        break;
      case TENSOR_EXCLD_2D_ONE:
				System.err.println("RSFExclusionZones: SINGLE TENSOR-ELLIPSE EXCLUDER");
				/////////////////////////////////////////////////
				// get parameters
				pfile = parms.getString("pfile","");
				if(pfile.equals("")) {
					System.err.println("A priority-map file was not specified.");
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
				c1 = parms.getFloat("c1",0.f);
				System.err.println("c1 = " + c1);
				c2 = parms.getFloat("c2",0.f);
				System.err.println("c2 = " + c2);


				/////////////////////////////////////////////////
				// read image file to mark with zones
				image2D = new float[n2][n1];
				System.err.println("Reading image to be marked.");
				input.read(image2D);



				/////////////////////////////////////////////////
				// read priority file
				pri2D = new float[n2][n1];
				Input spriInput2 = new Input("pfile");
				System.err.println("Reading priority-map file: " + pfile);
				spriInput2.read(pri2D);


				/////////////////////////////////////////////////
				// read tensors file
				float[][][] set = new float[4][n2][n1];
				Input stenInput = new Input("tfile");
				System.err.println("Reading tensor file: " + tfile);
				stenInput.read(set);


				/////////////////////////////////////////////////
				// create poi array
				poi = new PointOfInterest[1];
				indices = new int[2];
				indices[0] = (int)((c1-ox1)/dx1 + 0.5f);
				indices[1] = (int)((c2-ox2)/dx2 + 0.5f);
				System.err.println("pri2D.length    = " + pri2D.length);
				System.err.println("pri2D[0].length = " +pri2D[0].length);
				System.err.println("i2 = " + indices[1]);
				System.err.println("i1 = " + indices[0]);
				System.err.println("c2 = " + c2);
				System.err.println("c1 = " + c1);
				System.err.println("o2 = " + ox2);
				System.err.println("o1 = " + ox1);
				System.err.println("d2 = " + dx2);
				System.err.println("d1 = " + dx1);
				sval = pri2D[indices[1]][indices[0]];
				poi[0] = new PointOfInterest(sval,2,indices);


				/////////////////////////////////////////////////
				// create tensors
				System.err.println("Loading tensors.");
				EigenTensors2 set2 = new EigenTensors2(n1,n2);
				float[] su = new float[2];
				for(int i2=0;i2<n2;++i2)
					for(int i1=0;i1<n1;++i1) {
						su[0] = set[0][i2][i1];
						su[1] = set[1][i2][i1];
						set2.setEigenvectorU(i1,i2,su);
					}
				set2.setEigenvalues(set[2],set[3]);


				/////////////////////////////////////////////////////////
				// window picks and create image with exclusion zones
				System.err.println("Marking exclusion zone.");
				TensorEllipseCIPExcluder steExcluder2 = new TensorEllipseCIPExcluder(set2,rmax,n1,n2,false);
				zimage2D = ExclusionZones.showZones(image2D,poi,steExcluder2);


				/////////////////////////////////////////////////////////
				// Write image file and close inputs
				System.err.println("Writing image marked with exclusion zone.");
				output.setN(1,n1);
				output.setN(2,n2);
				output.setN(3,1);//rsf problem with 2D images keeping 3D data info
				output.setDelta(1,dx1);
				output.setDelta(2,dx2);
				output.setDelta(3,1.f);//rsf problem with 2D images keeping 3D data info
				output.setOrigin(1,ox1);
				output.setOrigin(2,ox2);
				output.setOrigin(3,0.f);//rsf problem with 2D images keeping 3D data info
				output.write(zimage2D);
				output.close();
				input.close();
				spriInput2.close();
				stenInput.close();
        break;
      case TENSOR_EXCLD_3D_ONE:
				System.err.println("RSFExclusionZones: SINGLE TENSOR-ELLIPSOID EXCLUDER");
				/////////////////////////////////////////////////
				// get parameters
				pfile = parms.getString("pfile","");
				if(pfile.equals("")) {
					System.err.println("A priority-map file was not specified.");
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
				c1 = parms.getFloat("c1",0.f);
				System.err.println("c1 = " + c1);
				c2 = parms.getFloat("c2",0.f);
				System.err.println("c2 = " + c2);
				c3 = parms.getFloat("c3",0.f);
				System.err.println("c3 = " + c3);


				/////////////////////////////////////////////////
				// read image file to mark with zones
				image3D = new float[n3][n2][n1];
				System.err.println("Reading image to be marked.");
				input.read(image3D);


				/////////////////////////////////////////////////
				// read priority file
				pri3D = new float[n3][n2][n1];
				Input stpriInput3D = new Input("pfile");
				System.err.println("Reading priority-map file: " + pfile);
				stpriInput3D.read(pri3D);


				/////////////////////////////////////////////////
				// read tensors file
				float[][][][] set3D = new float[9][n3][n2][n1];
				System.err.println("Reading tensor file: " + tfile);
				Input stenInput3D = new Input("tfile");
				stenInput3D.read(set3D);


				/////////////////////////////////////////////////
				// create poi array
				poi = new PointOfInterest[1];
				indices = new int[3];
				indices[0] = (int)((c1-ox1)/dx1 + 0.5f);
				indices[1] = (int)((c2-ox2)/dx2 + 0.5f);
				indices[2] = (int)((c3-ox3)/dx3 + 0.5f);
				sval = pri3D[indices[2]][indices[1]][indices[0]];
				poi[0] = new PointOfInterest(sval,3,indices);


				/////////////////////////////////////////////////
				// create tensors
				System.err.println("Loading tensors.");
				EigenTensors3 set3 = new EigenTensors3(n1,n2,n3,true);
				float su1,su2,su3,sw1,sw2,sw3;
				for(int i3=0;i3<n3;++i3)
					for(int i2=0;i2<n2;++i2)
						for(int i1=0;i1<n1;++i1) {
						  su1 = set3D[0][i3][i2][i1];
						  su2 = set3D[1][i3][i2][i1];
						  su3 = set3D[2][i3][i2][i1];
						  sw1 = set3D[3][i3][i2][i1];
						  sw2 = set3D[4][i3][i2][i1];
						  sw3 = set3D[5][i3][i2][i1];
						  set3.setEigenvectorU(i1,i2,i3,su1,su2,su3);
						  set3.setEigenvectorW(i1,i2,i3,sw1,sw2,sw3);
						}
				set3.setEigenvalues(set3D[6],set3D[7],set3D[8]);


				/////////////////////////////////////////////////////////
				// window picks and create image with exclusion zones
				System.err.println("Marking exclusion zone.");
				TensorEllipsoidCIPExcluder steExcluder3 = new TensorEllipsoidCIPExcluder(set3,rmax,n1,n2,n3,false);
				zimage3D = ExclusionZones.showZones(image3D,poi,steExcluder3);


				/////////////////////////////////////////////////////////
				// Write image file and close inputs
				System.err.println("Writing image with exclusion zone.");
				output.setN(1,n1);
				output.setN(2,n2);
				output.setN(3,n3);
				output.setDelta(1,dx1);
				output.setDelta(2,dx2);
				output.setDelta(3,dx3);
				output.setOrigin(1,ox1);
				output.setOrigin(2,ox2);
				output.setOrigin(3,ox3);
				output.write(zimage3D);
				output.close();
				input.close();
				stpriInput3D.close();
				stenInput3D.close();
        break;
      case ELLP_EXCLD_2D:
				System.err.println("RSFExclusionZones: ELLIPSE EXCLUDER");
				/////////////////////////////////////////////////
				// get parameters
				cfile = parms.getString("cfile","");
				if(cfile.equals("")) {
					System.err.println("A coordinate file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				pfile = parms.getString("pfile","");
				if(pfile.equals("")) {
					System.err.println("A priority-map file was not specified.");
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
				f = parms.getInt("first",0);
				j = parms.getInt("jump",1);


				/////////////////////////////////////////////////
				// read image file to mark with zones
				image2D = new float[n2][n1];
				System.err.println("Reading image to be marked.");
				input.read(image2D);

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
				// read priority file
				pri2D = new float[n2][n1];
				Input priInput = new Input("pfile");
				System.err.println("Reading priority-map file: " + pfile);
				priInput.read(pri2D);

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
				// window picks and create image with exclusion zones
				System.err.println("Marking exclusion zones.");
				EllipseCIPExcluder eExcluder2 = new EllipseCIPExcluder(r1,r2,n1,n2,false);
				wpoi = WindowCIPPicks.window(poi,f,poi.length-1,j,poi.length);
				System.err.println("Number of possible zones = " + poi.length);
				System.err.println("Number of zones marked = " + wpoi.length);
				zimage2D = ExclusionZones.showZones(image2D,wpoi,eExcluder2);


				/////////////////////////////////////////////////////////
				// Write image file and close inputs
				System.err.println("Writing image marked with exclusion zones.");
				output.setN(1,n1);
				output.setN(2,n2);
				output.setN(3,1);//rsf problem with 2D images keeping 3D data info
				output.setDelta(1,dx1);
				output.setDelta(2,dx2);
				output.setDelta(3,1.f);//rsf problem with 2D images keeping 3D data info
				output.setOrigin(1,ox1);
				output.setOrigin(2,ox2);
				output.setOrigin(3,0.f);//rsf problem with 2D images keeping 3D data info
				output.write(zimage2D);
				output.close();
				input.close();
				cInput.close();
				priInput.close();
        break;
      case ELLP_EXCLD_3D:
				System.err.println("RSFExclusionZones: ELLIPSOID EXCLUDER");
				/////////////////////////////////////////////////
				// get parameters
				cfile = parms.getString("cfile","");
				if(cfile.equals("")) {
					System.err.println("A coordinate file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				pfile = parms.getString("pfile","");
				if(pfile.equals("")) {
					System.err.println("A priority-map file was not specified.");
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
				f = parms.getInt("first",0);
				j = parms.getInt("jump",1);


				/////////////////////////////////////////////////
				// read image file to mark with zones
				image3D = new float[n3][n2][n1];
				System.err.println("Reading image to be marked.");
				input.read(image3D);

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
				// read priority file
				pri3D = new float[n3][n2][n1];
				Input priInput3D = new Input("pfile");
				System.err.println("Reading priority-map file: " + pfile);
				priInput3D.read(pri3D);

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
				// window picks and create image with exclusion zones
				System.err.println("Marking exclusion zones.");
				EllipsoidCIPExcluder eExcluder3 = new EllipsoidCIPExcluder(r1,r2,r3,n1,n2,n3,false);
				wpoi = WindowCIPPicks.window(poi,f,npicks-1,j,npicks);
				System.err.println("Number of possible zones = " + poi.length);
				System.err.println("Number of zones marked = " + wpoi.length);
				zimage3D = ExclusionZones.showZones(image3D,wpoi,eExcluder3);


				/////////////////////////////////////////////////////////
				// Write image file and close inputs
				System.err.println("Writing image marked with exclusion zones.");
				output.setN(1,n1);
				output.setN(2,n2);
				output.setN(3,n3);
				output.setDelta(1,dx1);
				output.setDelta(2,dx2);
				output.setDelta(3,dx3);
				output.setOrigin(1,ox1);
				output.setOrigin(2,ox2);
				output.setOrigin(3,ox3);
				output.write(zimage3D);
				output.close();
				input.close();
				cInput3D.close();
				priInput3D.close();
        break;
      case TENSOR_EXCLD_2D:
				System.err.println("RSFExclusionZones: TENSOR-ELLIPSE EXCLUDER");
				/////////////////////////////////////////////////
				// get parameters
				cfile = parms.getString("cfile","");
				if(cfile.equals("")) {
					System.err.println("A coordinate file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				pfile = parms.getString("pfile","");
				if(pfile.equals("")) {
					System.err.println("A priority-map file was not specified.");
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
				f = parms.getInt("first",0);
				j = parms.getInt("jump",1);


				/////////////////////////////////////////////////
				// read image file to mark with zones
				image2D = new float[n2][n1];
				System.err.println("Reading image to be marked.");
				input.read(image2D);


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
				Input priInput2 = new Input("pfile");
				System.err.println("Reading priority-map file: " + pfile);
				priInput2.read(pri2D);


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
				// window picks and create image with exclusion zones
				System.err.println("Marking exclusion zones.");
				System.err.println("rmax = " + rmax);
				TensorEllipseCIPExcluder teExcluder2 = new TensorEllipseCIPExcluder(et2,rmax,n1,n2,false);
				wpoi = WindowCIPPicks.window(poi,f,poi.length-1,j,poi.length);
				System.err.println("Number of possible zones = " + poi.length);
				System.err.println("Number of zones marked = " + wpoi.length);
				zimage2D = ExclusionZones.showZones(image2D,wpoi,teExcluder2);


				/////////////////////////////////////////////////////////
				// Write image file and close inputs
				System.err.println("Writing image marked with exclusion zones.");
				output.setN(1,n1);
				output.setN(2,n2);
				output.setN(3,1);//rsf problem with 2D images keeping 3D data info
				output.setDelta(1,dx1);
				output.setDelta(2,dx2);
				output.setDelta(3,1.f);//rsf problem with 2D images keeping 3D data info
				output.setOrigin(1,ox1);
				output.setOrigin(2,ox2);
				output.setOrigin(3,0.f);//rsf problem with 2D images keeping 3D data info
				output.write(zimage2D);
				output.close();
				input.close();
				cInput2.close();
				priInput2.close();
				tenInput.close();
        break;
      case TENSOR_EXCLD_3D:
				System.err.println("RSFExclusionZones: TENSOR-ELLIPSOID EXCLUDER");
				/////////////////////////////////////////////////
				// get parameters
				cfile = parms.getString("cfile","");
				if(cfile.equals("")) {
					System.err.println("A coordinate file was not specified.");
					output.close();
					input.close();
					System.exit(-1);
				}
				pfile = parms.getString("pfile","");
				if(pfile.equals("")) {
					System.err.println("A priority-map file was not specified.");
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
				f = parms.getInt("first",0);
				j = parms.getInt("jump",1);


				/////////////////////////////////////////////////
				// read image file to mark with zones
				image3D = new float[n3][n2][n1];
				System.err.println("Reading image to be marked.");
				input.read(image3D);


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
				// read priority file
				pri3D = new float[n3][n2][n1];
				Input tpriInput3D = new Input("pfile");
				System.err.println("Reading priority-map file: " + pfile);
				tpriInput3D.read(pri3D);


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
				// window picks and create image with exclusion zones
				System.err.println("Marking exclusion zones.");
				TensorEllipsoidCIPExcluder teExcluder3 = new TensorEllipsoidCIPExcluder(et3,rmax,n1,n2,n3,false);
				wpoi = WindowCIPPicks.window(poi,f,npicks-1,j,npicks);
				System.err.println("Number of possible zones = " + poi.length);
				System.err.println("Number of zones marked = " + wpoi.length);
				zimage3D = ExclusionZones.showZones(image3D,wpoi,teExcluder3);


				/////////////////////////////////////////////////////////
				// Write image file and close inputs
				System.err.println("Writing image with exclusion zones.");
				output.setN(1,n1);
				output.setN(2,n2);
				output.setN(3,n3);
				output.setDelta(1,dx1);
				output.setDelta(2,dx2);
				output.setDelta(3,dx3);
				output.setOrigin(1,ox1);
				output.setOrigin(2,ox2);
				output.setOrigin(3,ox3);
				output.write(zimage3D);
				output.close();
				input.close();
				tcInput3D.close();
				tpriInput3D.close();
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
			if(sopt_str.equals("tsingle"))
				return SwitchOpt.TENSOR_EXCLD_3D_ONE;
			if(sopt_str.equals("tensor"))
				return SwitchOpt.TENSOR_EXCLD_3D;
			if(sopt_str.equals("esingle"))
				return SwitchOpt.ELLP_EXCLD_3D_ONE;
			if(sopt_str.equals("ellp"))
				return SwitchOpt.ELLP_EXCLD_3D;
			return SwitchOpt.DEFAULT;
		}
		if(sopt_str.equals("tsingle"))
			return SwitchOpt.TENSOR_EXCLD_2D_ONE;
		if(sopt_str.equals("tensor"))
			return SwitchOpt.TENSOR_EXCLD_2D;
		if(sopt_str.equals("esingle"))
			return SwitchOpt.ELLP_EXCLD_2D_ONE;
		if(sopt_str.equals("ellp"))
			return SwitchOpt.ELLP_EXCLD_2D;
		return SwitchOpt.DEFAULT;
	} //end getSwitchOpt()

} //end class
