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
 * The class will compute the semblance within exclusion zones and return
 * a POI array sorted by the semblance value for each respective pick. 
 * The val for each POI returned is set to the exclusion-zone-semblance 
 * value.
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.11.14
 */
public class EZSemblanceSort
{

  /**
   * Zeros voxels within exclusion zones.
   * @param f a 2D priority map.
   * @param poi an array of PointOfInterest objects correspoinding to 
	 * CIP/image locations.
   * @param ezSemblance a CIPEZSemblance use to zero image voxels withing an
	 * exclusion zone.
   */
	public static void sortByEZSemblance(float[][] f, 
																			 PointOfInterest[] poi,
																			 CIPEZSemblance2 ezSemblance) { 

		int n1 = f[0].length;
		int n2 = f.length;
		PointOfInterest[][] poi2D = new PointOfInterest[n2][n1];
		initPOIArray(f,poi2D);
		ezSemblance.setPOIArray(poi2D);

		float val = 0.f;
		for(int i=0; i<poi.length; ++i) {
			val = ezSemblance.semblanceAt(poi[i].getIndex(0),
																	  poi[i].getIndex(1));
		  poi[i].setVal(val);
		}

		Arrays.sort(poi);
	}

  /**
   * Zeros voxels within exclusion zones.
   * @param f a 3D priority map.
   * @param poi an array of PointOfInterest objects correspoinding to 
	 * CIP/image locations.
   * @param ezSemblance a CIPEZSemblance use to zero image voxels withing an
	 * exclusion zone.
   */
	public static void sortByEZSemblance(float[][][] f, 
																			 PointOfInterest[] poi,
																			 CIPEZSemblance3 ezSemblance) { 

		int n1 = f[0][0].length;
		int n2 = f[0].length;
		int n3 = f.length;
		PointOfInterest[][][] poi3D = new PointOfInterest[n3][n2][n1];
		initPOIArray(f,poi3D);
		ezSemblance.setPOIArray(poi3D);

		float val = 0.f;
		for(int i=0; i<poi.length; ++i) {
			val = ezSemblance.semblanceAt(poi[i].getIndex(0),
																	  poi[i].getIndex(1),
																	  poi[i].getIndex(2));
		  poi[i].setVal(val);
		}

		Arrays.sort(poi);
	}


	////////////////////////////////////////////////////////////////////////////
	// Private

	private static void initPOIArrayLin2(PointOfInterest[] opoi, PointOfInterest[] ipoi) {
    int[] indices = new int[2];
		PointOfInterest tpoi;
    for(int i=0; i<ipoi.length; ++i) {
        indices[0] = ipoi[i].getIndex(0); 
        indices[1] = ipoi[i].getIndex(1); 
		    float val = ipoi[i].getVal();
				opoi[i] = new PointOfInterest(val,2,indices);
		}
	}

	private static void initPOIArrayLin3(PointOfInterest[] opoi, PointOfInterest[] ipoi) {
    int[] indices = new int[3];
		PointOfInterest tpoi;
    for(int i=0; i<ipoi.length; ++i) {
        indices[0] = ipoi[i].getIndex(0); 
        indices[1] = ipoi[i].getIndex(1); 
        indices[2] = ipoi[i].getIndex(2); 
		    float val = ipoi[i].getVal();
				opoi[i] = new PointOfInterest(val,3,indices);
		}
	}

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
