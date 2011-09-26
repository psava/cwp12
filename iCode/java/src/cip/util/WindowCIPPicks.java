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
 * Provides functions for windowing CIP picks/locations.
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.11.14
 */
public class WindowCIPPicks
{

  /**
   * Windows CIP picks.
   * @param ipoi an array of PointOfInterest objects correspoinding to 
	 * CIP picks/locations.
   * @param f the first CIP pick/location to start from.
   * @param l the last CIP pick/location to pick from.
   * @param j jump j locations between picks/locations.
   * @param n the maximum number CIP picks/locations.
   * @return PointOfInterest[] of windowed picks/locations.
   */
	public static PointOfInterest[] window(PointOfInterest[] ipoi, 
																				 int f, int l, int j, int n) {
	  
		assert 0 <= f;
		assert n <= ipoi.length;
		assert l < n;
		assert 1 <= j;
		int npoi = (l - f + 1)/j;
		if((l - f + 1)%j > 0)
		  ++npoi;

		PointOfInterest[] opoi;
		if(npoi <= n)
			opoi = new PointOfInterest[npoi];
		else
			opoi = new PointOfInterest[n];

		int io = 0;
		for(int ii=f; ii<=l && io<n; ii+=j,++io)
		  opoi[io] = ipoi[ii];

	  return opoi;
	}
  /**

   * Windows CIP picks at random.
   * @param ipoi an array of PointOfInterest objects correspoinding to 
	 * CIP picks/locations.
   * @param f the first CIP pick/location to start from.
   * @param l the last CIP pick/location to pick from.
   * @param j jump j locations between picks/locations.
   * @return PointOfInterest[] of windowed picks/locations.
   */
	public static PointOfInterest[] windowRandom(PointOfInterest[] ipoi, 
																							 int f, int l, int n) {
	  
	  return ipoi;
	}

}
