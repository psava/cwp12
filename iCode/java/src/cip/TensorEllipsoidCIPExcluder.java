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
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.11.14
 */
public class TensorEllipsoidCIPExcluder implements CIPExcluder3
{

  /**
   * Constructs a Tensor based non-euclidean ellipsoid CIP excluder.  NOTE:
	 * When using this constructor the setPOIArray() function must be called
	 * before any other function that opperations on this class.
	 * @param et an EigenTensors3 field for a 3D image.
   * @param maxr the maximum radius distance for the exclusion zone.
   * @param n1 fastest dim size.
   * @param n2 the next fastest dim size.
   * @param n3 slowest dim size.
   * @param threaded a boolean. If true: locations are excluded in parallel.
   */
	public TensorEllipsoidCIPExcluder(EigenTensors3 et, float maxr, 
																	  int n1, int n2, int n3,
																	  boolean threaded) {
	  _maxr = maxr;
		_n1 = n1;
		_n2 = n2;
		_n3 = n3;
		_et = et;
		_threaded = threaded;
		if(!_threaded) {
			_ete = new ExclusionTensorEllipsoid(_n1,_n2,_n3,_et);
			_ete.setSize(_maxr);
		}
		_poiArray = null;
	}

  /**
   * Constructs a Tensor based non-euclidean ellipsoid CIP excluder.
	 * @param et an EigenTensors3 field for a 3D image.
   * @param maxr the maximum radius distance for the exclusion zone.
   * @param n1 fastest dim size.
   * @param n2 the next fastest dim size.
   * @param n3 slowest dim size.
   * @param poiArray A 3D PointOfInterest array the corresponds to a 3D
	 * priority map.
   * @param threaded a boolean. If true: locations are excluded in parallel.
   */
	public TensorEllipsoidCIPExcluder(EigenTensors3 et, float maxr,
																	  int n1, int n2, int n3,
																	  PointOfInterest[][][] poiArray, 
																	  boolean threaded) {
	  _maxr = maxr;
		_n1 = n1;
		_n2 = n2;
		_n3 = n3;
		_et = et;
		_threaded = threaded;
		if(!_threaded) {
			_ete = new ExclusionTensorEllipsoid(_n1,_n2,_n3,_et);
			_ete.setSize(_maxr);
		}
		_poiArray = poiArray;
	}

  /**
   * Excludes points around the index i1,i2,i3.
   * @param i1 index in the fastest dim.
   * @param i2 index in the next fastest dim.
   * @param i3 index in the slowest dim.
   */
	public void excludeAt(int i1, int i2, int i3) {
		if(_poiArray == null) {
		  System.err.println("A PointOfInterest Array was not been declared.");
			System.exit(-1);
		}
	  if(_threaded)
			excludeThreaded(i1, i2, i3);
		else
			excludeLoop(i1, i2, i3);
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
   * Excludes points around the index i1,i2,i3 using a serial loop.
   */
	private void excludeLoop(int i1, int i2, int i3) {
		_ete.setLocation(i1,i2,i3);
    PointOfInterest tpoi;
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
						tpoi.setVal(0.f);
						tpoi.setLock(true);
						tpoi.setLive(false);
					}
				}
			}
		}
	}

  /**
   * Excludes points around the index i1,i2,i3 using multible threads.
   */
	private void excludeThreaded(final int i1, final int i2, final int i3) {
		int _r = (int) (_maxr + 0.5f);
		final int bx1 = i1-_r > 0 ? i1-_r : 0;
		final int bx2 = i2-_r > 0 ? i2-_r : 0;
		final int bx3 = i3-_r > 0 ? i3-_r : 0;
		final int ex1 = i1+_r < _n1 ? i1+_r : _n1-1;
		final int ex2 = i2+_r < _n2 ? i2+_r : _n2-1;
		final int ex3 = i3+_r < _n3 ? i3+_r : _n3-1;
    final AtomicInteger ai = new AtomicInteger(bx3);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
					PointOfInterest tpoi;
					float d;
					ExclusionTensorEllipsoid _eteThread = new ExclusionTensorEllipsoid(_n1,_n2,_n3,_et);
					_eteThread.setSize(_maxr);
					_eteThread.setLocation(i1,i2,i3);
          for (int ii3=ai.getAndIncrement(); ii3<=ex3; ii3=ai.getAndIncrement()) {
						for(int ii2=bx2; ii2<=ex2; ++ii2) {
							for(int ii1=bx1; ii1<=ex1; ++ii1) {
								tpoi = _poiArray[ii3][ii2][ii1];
									if(tpoi.isLocked())
										continue;
								d = _eteThread.getDistance(ii1,ii2,ii3);
								if(d <= _maxr) {
									tpoi.setVal(0.f);
									tpoi.setLock(true);
									tpoi.setLive(false);
								}
							}
						}
					}
				}
			});
		}
    Threads.startAndJoin(threads);
	}

	private int _n1,_n2,_n3;
	private float _maxr;
	boolean _threaded;
	private EigenTensors3 _et;
	private ExclusionTensorEllipsoid _ete;
  private PointOfInterest[][][] _poiArray;
}
