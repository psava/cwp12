/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

//  GreedyCIPPicker2.java
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

public class GreedyCIPPicker2
{

	public GreedyCIPPicker2(float[][] pMap, CIPExcluder2 excluder) {
    _n2 = pMap.length;
    _n1 = pMap[0].length;
    _poiArray = new PointOfInterest[_n2][_n1];
		_excluder = excluder;
    _poiPickFrom = new PointOfInterest[_n2*_n1];
	  _poiPQPicks = new PriorityBlockingQueue<PointOfInterest>();
    initPOIArrayAndQueue(pMap);
		_excluder.setPOIArray(_poiArray);
		Arrays.sort(_poiPickFrom);
    _numpicks = pickCIPs();
    _poiPicks = new PointOfInterest[1];
    _poiPicks = _poiPQPicks.toArray(_poiPicks);
		Arrays.sort(_poiPicks);
	}

	public int[][] getCIPIndices() {
    int[][] indices = new int[2][_numpicks];
    for(int i=0; i<_numpicks; ++i) {
      indices[0][i] = _poiPicks[i].getIndex(0);
      indices[1][i] = _poiPicks[i].getIndex(1);
    }
    return indices;
	}

	public float[][] getCIPCoord(float ox1, float ox2, 
															 float dx1, float dx2) {
    float[][] coord = new float[2][_numpicks];
    for(int i=0; i<_numpicks; ++i) {
      coord[0][i] = ox1 + dx1*_poiPicks[i].getIndex(0);
      coord[1][i] = ox2 + dx2*_poiPicks[i].getIndex(1);
    }
    return coord;
	}

	public PointOfInterest[] getLinearPOIPicksArray() { return _poiPicks; }

	public PointOfInterest[][] getPOIArray() { return _poiArray; }

	public int getNumPicks() { return _numpicks; }



/////////////////////////////////////////////////////////////
// Private

	private int pickCIPs() {
		int[] loc = new int[2];
	  PointOfInterest tpoi;
		for(int i=0; i<_n1*_n2; ++i) {
			tpoi = _poiPickFrom[i];
			if(tpoi.isLocked())
			  continue;
			if(tpoi.getVal() <= 0) {
				tpoi.setLock(true);
				tpoi.setLive(false);
			  continue;
			}
			_poiPQPicks.add(tpoi);
			tpoi.setLock(true);
			int i1 = tpoi.getIndex(0);
			int i2 = tpoi.getIndex(1);
			_excluder.excludeAt(i1,i2);
		}
    return _poiPQPicks.size();
	}

	private void initPOIArrayAndQueue(float[][] pMap) {
    int[] indices = new int[2];
    for(int i2=0; i2<_n2; ++i2) {
			int stride = _n1*i2;
      for(int i1=0; i1<_n1; ++i1) {
        indices[0] = i1;
        indices[1] = i2;
		    float val = pMap[i2][i1];
		    PointOfInterest tpoi = new PointOfInterest(val,2,indices);
		    _poiArray[i2][i1] = tpoi;
				_poiPickFrom[i1+stride] = tpoi;
      }
    }
	}

	private int _n1,_n2;
	private int _numpicks;
	private CIPExcluder2 _excluder;
  private PointOfInterest[][] _poiArray;
  private PointOfInterest[] _poiPicks;
  private PointOfInterest[] _poiPickFrom;
	private PriorityBlockingQueue<PointOfInterest> _poiPQPicks;
}
