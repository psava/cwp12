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
 * Uses TimeSolver2 (Ikonal Solver) to determine the exclusion zones.
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.11.14
 */
public class TensorEllipseCIPExcluder implements CIPExcluder2
{

  /**
   * Constructs a Tensor based non-euclidean ellipse CIP excluder.  NOTE:
	 * When using this constructor the setPOIArray() function must be called
	 * before any other function that opperations on this class.
	 * @param et an EigenTensors2 field for a 2D image.
   * @param maxr the maximum radius distance for the exclusion zone.
   * @param n1 fastest dim size.
   * @param n2 slowest dim size.
   * @param threaded a boolean. If true: locations are excluded in parallel.
   */
	public TensorEllipseCIPExcluder(EigenTensors2 et, float maxr, 
																  int n1, int n2,
																  boolean threaded) {
	  _maxr = maxr;
		_n1 = n1;
		_n2 = n2;
		_et = et;
		_threaded = threaded;
		if(!_threaded) {
			_ete = new ExclusionTensorEllipse(_n1,_n2,_et);
			_ete.setSize(_maxr);
		}
		_poiArray = null;
	}

  /**
   * Constructs a Tensor based non-euclidean ellipse CIP excluder.
	 * @param et an EigenTensors2 field for a 2D image.
   * @param maxr the maximum radius distance for the exclusion zone.
   * @param n1 fastest dim size.
   * @param n2 slowest dim size.
   * @param poiArray A 2D PointOfInterest array the corresponds to a 2D
	 * priority map.
   * @param threaded a boolean. If true: locations are excluded in parallel.
   */
	public TensorEllipseCIPExcluder(EigenTensors2 et, float maxr,
																  int n1, int n2,
																  PointOfInterest[][] poiArray, 
																  boolean threaded) {
	  _maxr = maxr;
		_n1 = n1;
		_n2 = n2;
		_et = et;
		_threaded = threaded;
		if(!_threaded) {
			_ete = new ExclusionTensorEllipse(_n1,_n2,_et);
			_ete.setSize(_maxr);
		}
		_poiArray = poiArray;
	}

  /**
   * Excludes points around the index i1,i2.
   * @param i1 index in the fastest dim.
   * @param i2 index in the slowest dim.
   */
	public void excludeAt(int i1, int i2) {
		if(_poiArray == null) {
		  System.err.println("A PointOfInterest Array was not been declared.");
			System.exit(-1);
		}
	  if(_threaded)
			excludeThreaded(i1,i2);
		else
			excludeLoop(i1,i2);
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
   * Excludes points around the index i1,i2 using a serial loop.
   */
	private void excludeLoop(int i1, int i2) {
		_ete.setLocation(i1,i2);
    PointOfInterest tpoi;
    float d;
		int _r = (int) (_maxr + 0.5f);
		int bx1 = i1-_r > 0 ? i1-_r : 0;
		int bx2 = i2-_r > 0 ? i2-_r : 0;
		int ex1 = i1+_r < _n1 ? i1+_r : _n1-1;
		int ex2 = i2+_r < _n2 ? i2+_r : _n2-1;
		for(int ii2=bx2; ii2<=ex2; ++ii2) {
			for(int ii1=bx1; ii1<=ex1; ++ii1) {
				tpoi = _poiArray[ii2][ii1];
				if(tpoi.isLocked())
					continue;
				d = _ete.getDistance(ii1,ii2);
				if(d <= _maxr) {
					tpoi.setVal(0.f);
					tpoi.setLock(true);
					tpoi.setLive(false);
				}
			}
		}
	}

  /**
   * Excludes points around the index i1,i2 using multible threads.
   */
	private void excludeThreaded(final int i1, final int i2) {
		int _r = (int) (_maxr + 0.5f);
		final int bx1 = i1-_r > 0 ? i1-_r : 0;
		final int bx2 = i2-_r > 0 ? i2-_r : 0;
		final int ex1 = i1+_r < _n1 ? i1+_r : _n1-1;
		final int ex2 = i2+_r < _n2 ? i2+_r : _n2-1;
    final AtomicInteger ai = new AtomicInteger(bx2);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
					PointOfInterest tpoi;
					float d;
					ExclusionTensorEllipse _eteThread = new ExclusionTensorEllipse(_n1,_n2,_et);
					_eteThread.setSize(_maxr);
					_eteThread.setLocation(i1,i2);
          for (int ii2=ai.getAndIncrement(); ii2<=ex2; ii2=ai.getAndIncrement()) {
						for(int ii1=bx1; ii1<=ex1; ++ii1) {
							tpoi = _poiArray[ii2][ii1];
								if(tpoi.isLocked())
									continue;
							d = _eteThread.getDistance(ii1,ii2);
							if(d <= _maxr) {
								tpoi.setVal(0.f);
								tpoi.setLock(true);
								tpoi.setLive(false);
							}
						}
					}
				}
			});
		}
    Threads.startAndJoin(threads);
	}

	private int _n1,_n2;
	private float _maxr;
	boolean _threaded;
	private EigenTensors2 _et;
	private ExclusionTensorEllipse _ete;
  private PointOfInterest[][] _poiArray;
}
