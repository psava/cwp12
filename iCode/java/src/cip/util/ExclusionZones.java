/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package cip.util;

import java.util.*;
import static java.lang.Math.sqrt;
import java.io.IOException;
import edu.mines.jtk.dsp.*;
import cip.*;

/**
 * Provided Factor Functions for Displaying exclusion zones on images.
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.11.14
 */
public class ExclusionZones
{

  /**
   * Zeros pixels within exclusion zones.
   * @param f an image of pixels.
   * @param poi an array of PointOfInterest objects correspoinding to 
	 * CIP/image locations.
   * @param excluder a CIPExcluder use to zero image pixels withing an
	 * exclusion zone.
   * @return float[][] an image of f where the pixels for each exclusion 
	 * zone at each location in poi are zeroed. 
   */
	public static float[][] showZones(float[][] f, 
																		PointOfInterest[] poi,
																		CIPExcluder2 excluder) { 

		int n1 = f[0].length;
		int n2 = f.length;
		float[][] of = new float[n2][n1];
		PointOfInterest[][] poi2D = new PointOfInterest[n2][n1];
		initPOIArray(f,poi2D);
		excluder.setPOIArray(poi2D);

		for(int i=0; i<poi.length; ++i) {
		  excluder.excludeAt(poi[i].getIndex(0),poi[i].getIndex(1));
		}

		for(int i2=0; i2<n2; ++i2)
			for(int i1=0; i1<n1; ++i1)
			  if(poi2D[i2][i1].isLocked())
				  of[i2][i1] = 0.f;
				else {
				  of[i2][i1] = poi2D[i2][i1].getVal();
				}
	
	  return of;
	}

  /**
   * Zeros voxels within exclusion zones.
   * @param f an image of voxels.
   * @param poi an array of PointOfInterest objects correspoinding to 
	 * CIP/image locations.
   * @param excluder a CIPExcluder use to zero image voxels withing an
	 * exclusion zone.
   * @return float[][][] an image of f where the voxels for each exclusion 
	 * zone at each location in poi are zeroed. 
   */
	public static float[][][] showZones(float[][][] f, 
																			PointOfInterest[] poi,
																			CIPExcluder3 excluder) { 

		int n1 = f[0][0].length;
		int n2 = f[0].length;
		int n3 = f.length;
		float[][][] of = new float[n3][n2][n1];
		PointOfInterest[][][] poi3D = new PointOfInterest[n3][n2][n1];
		initPOIArray(f,poi3D);
		excluder.setPOIArray(poi3D);

		for(int i=0; i<poi.length; ++i)
		  excluder.excludeAt(poi[i].getIndex(0),
												 poi[i].getIndex(1),
												 poi[i].getIndex(2));

		for(int i3=0; i3<n3; ++i3)
			for(int i2=0; i2<n2; ++i2)
				for(int i1=0; i1<n1; ++i1)
					if(poi3D[i3][i2][i1].isLocked()) {
						of[i3][i2][i1] = 0.f;
					}
					else {
						of[i3][i2][i1] = poi3D[i3][i2][i1].getVal();
					}
	
	  return of;
	}

  /**
   * Zeros voxels within exclusion zones.
   * @param poi an array of PointOfInterest objects correspoinding to 
	 * CIP/image locations.
   * @param poi2D a 2D array of PointOfInterest objects correspoinding to
	 * image pixels.
   * @param ee an ExclusionEllipse distance map used to compute exclusion zone.
   * @return float[][][] an image of f where the voxels for each exclusion 
	 * zone at each location in poi are zeroed. 
   */
	public static int[][] getCIPNearestMaxIndices2D(PointOfInterest[] poi,
																									PointOfInterest[][] poi2D,
																								  ExclusionEllipse ee,
																									int r1, int r2) {
		int n1 = poi2D[0].length;
		int n2 = poi2D.length;
    int[][] indices = new int[2][poi.length];
    for(int i=0; i<poi.length; ++i) {
		  int i1 = poi[i].getIndex(0);
		  int i2 = poi[i].getIndex(1);
			int mi1 = i1;
			int mi2 = i2;
			ee.setLocation(i1,i2);
			float mval = poi[i].getVal();
			float d1,d2;
			float d,dmin,val;
			dmin = 2.f;
			int bx1 = i1-r1 > 0 ? i1-r1 : 0;
			int bx2 = i2-r2 > 0 ? i2-r2 : 0;
			int ex1 = i1+r1 < n1 ? i1+r1 : n1-1;
			int ex2 = i2+r2 < n2 ? i2+r2 : n2-1;
			for(int ii2=bx2; ii2<=ex2; ++ii2) {
				for(int ii1=bx1; ii1<=ex1; ++ii1) {
					d = ee.getDistance(ii1,ii2);
					val = poi2D[ii2][ii1].getVal();
					if(val >= mval && d <= 1.f) {
						if(d < dmin) {
							mval = val;
						  dmin = d;
							mi1=ii1;
							mi2=ii2;
						}
					}
				}
			}
      indices[0][i] = mi1;
      indices[1][i] = mi2;
    }
    return indices;
	}


	////////////////////////////////////////////////////////////////////////////
	// Private

	private static void initPOIArray(float[][] f, PointOfInterest[][] poi2D) {
		int n1 = poi2D[0].length;
		int n2 = poi2D.length;
    int[] indices = new int[2];
    for(int i2=0; i2<n2; ++i2)
      for(int i1=0; i1<n1; ++i1) {
        indices[0] = i1;
        indices[1] = i2;
		    float val = f[i2][i1];
		    poi2D[i2][i1] = new PointOfInterest(val,2,indices);
      }
	}

	private static void initPOIArray(float[][][] f, PointOfInterest[][][] poi3D) {
		int n1 = poi3D[0][0].length;
		int n2 = poi3D[0].length;
		int n3 = poi3D.length;
    int[] indices = new int[3];
    for(int i3=0; i3<n3; ++i3)
    for(int i2=0; i2<n2; ++i2)
      for(int i1=0; i1<n1; ++i1) {
        indices[0] = i1;
        indices[1] = i2;
        indices[2] = i3;
		    float val = f[i3][i2][i1];
		    poi3D[i3][i2][i1] = new PointOfInterest(val,3,indices);
      }
	}

}
