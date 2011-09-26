/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package cip;

/**
 * Excludes locations based on grid aligend ellipses.
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.10.29
 */
public class ExclusionEllipse implements CIPDistanceMap2
{

  /**
   * Constructs an exclusion ellise for specified tensors.
   * @param r1 exlusion radius in the 1st dimension.
   * @param r2 exlusion radius in the 2st dimension.
   */
  public ExclusionEllipse(int r1, int r2) {
	  _r1 = r1;
		_r2 = r2;
		_d1 = 0.f;
		_d2 = 0.f;
		_pii2 = -1;
  }

  /**
   * Sets the location (i1,i2) of this brush.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   */
  public void setLocation(int i1, int i2) {
	  _i1=i1;
	  _i2=i2;
  }


  /**
   * Gets the distance d that the point (ii1,ii2) is 
	 * from the current location.
   * @param ii1 index in 1st dimension.
   * @param ii2 index in 2nd dimension.
   */
  public float getDistance(int ii1, int ii2) {
		if(ii2 != _pii2) {
		  _pii2 = ii2;
			_d2 = _i2-ii2;
			_d2 *= _d2;
			_d2 /= _r2*_r2;
		}
		_d1 = _i1-ii1;
		_d1 *= _d1;
		_d1 /= _r1*_r1;
		return _d2 + _d1;
  }

  private int _r1,_r2;
	private float _d2,_d1;
	private int _pii2;
  private int _i1,_i2;
}
