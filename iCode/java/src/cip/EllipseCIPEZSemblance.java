/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

//
//  THOMAS CULLISON
//  Colorado School of Mines
//  Center for Wave Phenomena

package cip;

import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import static java.lang.Math.sqrt;
import java.io.IOException;

import edu.mines.jtk.util.Threads;

public class EllipseCIPEZSemblance implements CIPEZSemblance2
{

  /**
   * Constructs a grid aligned ellipse CIP exclusion
	 * zone semblance calculator. Note: When using this constructor 
	 * the setPOIArray() function must be called before any other 
	 * function that opperations on this class.
	 * @param et an EigenTensors2 field for a 2D image.
   * @param r1 radius of the ellipse in the fastest dim.
   * @param r2 radius of the ellipse in the slowest dim.
   * @param n1 fastest dim size.
   * @param n2 the next fastest dim size.
   */
	public EllipseCIPEZSemblance(int r1, int r2,
															 int n1, int n2) {
	  _r1 = r1;
	  _r2 = r2;
	  _n1 = n1;
	  _n2 = n2;
    _ee = new ExclusionEllipse(_r1,_r2);
		_poiArray = null;
	}

  /**
   * Constructs a Tensor based non-euclidean ellipse CIP exclusion 
	 * zone semblance calculator.
	 * @param et an EigenTensors2 field for a 2D image.
   * @param r1 radius of the ellipse in the fastest dim.
   * @param r2 radius of the ellipse in the slowest dim.
   * @param n1 fastest dim size.
   * @param n2 the next fastest dim size.
   * @param poiArray A 2D PointOfInterest array the corresponds to a 2D
	 * priority map.
   */
	public EllipseCIPEZSemblance(int r1, int r2,
															 int n1, int n2,
															 PointOfInterest[][] poiArray) {
	  _r1 = r1;
	  _r2 = r2;
	  _n1 = n1;
	  _n2 = n2;
    _ee = new ExclusionEllipse(_r1,_r2);
		_poiArray = poiArray;
	}

  /**
   * Computes the semblance within the exclusion zone at the index i1,i2.
   * @param i1 index in the fastest dim.
   * @param i2 index in the slowest dim.
   * @return float: the semblance value for the exclusion zone.
   */
	public float semblanceAt(int i1, int i2) {
		if(_poiArray == null) {
		  System.err.println("A PointOfInterest Array was not been declared.");
			System.exit(-1);
		}
		return semblanceLoop(i1, i2);
	}

  /**
   * Sets the 2D PointOfInterest array.
   * @param poiArray A 2D PointOfInterest array the corresponds to a 2D
	 * priority map.
   */
	public void setPOIArray(PointOfInterest[][] poiArray) {
	  _poiArray = poiArray;
	}

	///////////////////////////////////////////////////////////////
	// Private 

  /**
   * Computes the semblance within the exclusion zone around 
	 * the index i1,i2 using a serial loop.
   */
	private float semblanceLoop(int i1, int i2) {
		_ee.setLocation(i1,i2);
    PointOfInterest tpoi;
		int count = 0;
		float num = 0.f;
		float den = 0.f;
    float d;
		int bx1 = i1-_r1 > 0 ? i1-_r1 : 0;
		int bx2 = i2-_r2 > 0 ? i2-_r2 : 0;
		int ex1 = i1+_r1 < _n1 ? i1+_r1 : _n1-1;
		int ex2 = i2+_r2 < _n2 ? i2+_r2 : _n2-1;
		for(int ii2=bx2; ii2<=ex2; ++ii2) {
			for(int ii1=bx1; ii1<=ex1; ++ii1) {
				tpoi = _poiArray[ii2][ii1];
				if(tpoi.isLocked())
					continue;
				d = _ee.getDistance(ii1,ii2);
				if(d <= 1.f) {
					++count;
					float val = tpoi.getVal();
					den += val*val; 
					num += val;
					tpoi.setLock(true);
					tpoi.setLive(false);
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

	private int _n1,_n2;
	private int _r1,_r2;
	private ExclusionEllipse _ee;
  private PointOfInterest[][] _poiArray;
}
