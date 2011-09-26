/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

//  EllipseCIPExcluder.java
//
//  THOMAS CULLISON
//  Colorado School of Mines
//  Center for Wave Phenomena

package cip;

import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.io.IOException;

import edu.mines.jtk.util.Threads;

public class EllipseCIPExcluder implements CIPExcluder2
{

	public EllipseCIPExcluder(int r1, int r2, 
													  int n1, int n2, 
													  boolean threaded) {
		_r1 = r1;
		_r2 = r2;
		_n1 = n1;
		_n2 = n2;
		_threaded = threaded;
    _ee = new ExclusionEllipse(_r1,_r2);
		_poiArray = null;
	}

	public EllipseCIPExcluder(int r1, int r2,
													  int n1, int n2, 
													  PointOfInterest[][] poiArray, 
													  boolean threaded) {
		_r1 = r1;
		_r2 = r2;
		_n1 = n1;
		_n2 = n2;
		_threaded = threaded;
    _ee = new ExclusionEllipse(_r1,_r2);
		_poiArray = poiArray;
	}

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

	public void setPOIArray(PointOfInterest[][] poiArray) {
	  _poiArray = poiArray;
	}

	///////////////////////////////////////////////////////////////
	// Private 

	private void excludeLoop(int i1, int i2) {
		_ee.setLocation(i1,i2);
    PointOfInterest tpoi;
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
					tpoi.setLock(true);
					tpoi.setLive(false);
				}
			}
		}
	}

	//FIXME: results in picks to close to each other
	private void excludeThreaded(final int i1, final int i2) {
		final int bx1 = i1-_r1 > 0 ? i1-_r1 : 0;
		final int bx2 = i2-_r2 > 0 ? i2-_r2 : 0;
		final int ex1 = i1+_r1 < _n1 ? i1+_r1 : _n1-1;
		final int ex2 = i2+_r2 < _n2 ? i2+_r2 : _n2-1;
    final AtomicInteger ai = new AtomicInteger(bx2);
    Thread[] threads = Threads.makeArray();
    for (int ithread=0; ithread<threads.length; ++ithread) {
      threads[ithread] = new Thread(new Runnable() {
        public void run() {
					PointOfInterest tpoi;
					ExclusionEllipse _eeThread = new ExclusionEllipse(_r1,_r2);
					_eeThread.setLocation(i1,i2);
					float d;
          for (int ii2=ai.getAndIncrement(); ii2<=ex2; ii2=ai.getAndIncrement()) {
						for(int ii1=bx1; ii1<=ex1; ++ii1) {
							tpoi = _poiArray[ii2][ii1];
								if(tpoi.isLocked())
									continue;
							d = _eeThread.getDistance(ii1,ii2);
							if(d <= 1.f) {
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
	private int _r1, _r2;
	private boolean _threaded;
	private ExclusionEllipse _ee;
  private PointOfInterest[][] _poiArray;
}
