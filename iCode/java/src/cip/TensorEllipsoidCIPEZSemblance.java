/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

package cip;

import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import static java.lang.Math.sqrt;
import java.io.IOException;

import edu.mines.jtk.util.Threads;
import edu.mines.jtk.dsp.*;

/**
 * Uses TimeSolver3 (Ikonal Solver) to determine the exclusion zones.
 * The semblance within and exclusion zone can be calculated.
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.11.14
 */
public class TensorEllipsoidCIPEZSemblance implements CIPEZSemblance3
{

  /**
   * Constructs a Tensor based non-euclidean ellipsoid CIP exclusion
	 * zone semblance calculator. Note: When using this constructor 
	 * the setPOIArray() function must be called before any other 
	 * function that opperations on this class.
	 * @param et an EigenTensors3 field for a 3D image.
   * @param maxr the maximum radius distance for the exclusion zone.
   * @param n1 fastest dim size.
   * @param n2 the next fastest dim size.
   * @param n3 slowest dim size.
   */
	public TensorEllipsoidCIPEZSemblance(EigenTensors3 et, float maxr, 
																	  int n1, int n2, int n3) {
	  _maxr = maxr;
		_n1 = n1;
		_n2 = n2;
		_n3 = n3;
		_et = et;
		_ete = new ExclusionTensorEllipsoid(_n1,_n2,_n3,_et);
		_ete.setSize(_maxr);
		_poiArray = null;
	}

  /**
   * Constructs a Tensor based non-euclidean ellipsoid CIP exclusion 
	 * zone semblance calculator.
	 * @param et an EigenTensors3 field for a 3D image.
   * @param maxr the maximum radius distance for the exclusion zone.
   * @param n1 fastest dim size.
   * @param n2 the next fastest dim size.
   * @param n3 slowest dim size.
   * @param poiArray A 3D PointOfInterest array the corresponds to a 3D
	 * priority map.
   */
	public TensorEllipsoidCIPEZSemblance(EigenTensors3 et, float maxr,
																	  int n1, int n2, int n3,
																	  PointOfInterest[][][] poiArray) {
	  _maxr = maxr;
		_n1 = n1;
		_n2 = n2;
		_n3 = n3;
		_et = et;
		_ete = new ExclusionTensorEllipsoid(_n1,_n2,_n3,_et);
		_ete.setSize(_maxr);
		_poiArray = poiArray;
	}

  /**
   * Computes the semblance within the exclusion zone at the index i1,i2,i3.
   * @param i1 index in the fastest dim.
   * @param i2 index in the next fastest dim.
   * @param i3 index in the slowest dim.
   * @return float: the semblance value for the exclusion zone.
   */
	public float semblanceAt(int i1, int i2, int i3) {
		if(_poiArray == null) {
		  System.err.println("A PointOfInterest Array was not been declared.");
			System.exit(-1);
		}
		return semblanceLoop(i1, i2, i3);
	}

  /**
   * Sets the 3D PointOfInterest array.
   * @param poiArray A 3D PointOfInterest array the corresponds to a 3D
	 * priority map.
   */
	public void setPOIArray(PointOfInterest[][][] poiArray) {
	  _poiArray = poiArray;
	}

	///////////////////////////////////////////////////////////////
	// Private 

  /**
   * Computes the semblance within the exclusion zone around 
	 * the index i1,i2,i3 using a serial loop.
   */
	private float semblanceLoop(int i1, int i2, int i3) {
		_ete.setLocation(i1,i2,i3);
    PointOfInterest tpoi;
		int count = 0;
		float num = 0.f;
		float den = 0.f;
    float d;
		int _r = (int) (_maxr + 0.5f);
		int bx1 = i1-_r > 0 ? i1-_r : 0;
		int bx2 = i2-_r > 0 ? i2-_r : 0;
		int bx3 = i3-_r > 0 ? i3-_r : 0;
		int ex1 = i1+_r < _n1 ? i1+_r : _n1-1;
		int ex2 = i2+_r < _n2 ? i2+_r : _n2-1;
		int ex3 = i3+_r < _n3 ? i3+_r : _n3-1;
		for(int ii3=bx3; ii3<=ex3; ++ii3) {
			for(int ii2=bx2; ii2<=ex2; ++ii2) {
				for(int ii1=bx1; ii1<=ex1; ++ii1) {
					tpoi = _poiArray[ii3][ii2][ii1];
					if(tpoi.isLocked())
					  continue;
					d = _ete.getDistance(ii1,ii2,ii3);
					if(d <= _maxr) {
						++count;
						float val = tpoi.getVal();
						den += val*val; 
						num += val;
						tpoi.setLock(true);
						tpoi.setLive(false);
					}
				}
			}
		}
		float scale = 1.f/count;
		num *= num;
		num *= scale;
		float rval = 0.f;
		if(num != 0 && den != 0)
			rval = num/den;
		if(rval > 1.f) // roundoff can least to 1 + esp vals
		  rval = 1.f;
		return rval;
	}

	private int _n1,_n2,_n3;
	private float _maxr;
	private EigenTensors3 _et;
	private ExclusionTensorEllipsoid _ete;
  private PointOfInterest[][][] _poiArray;
}
