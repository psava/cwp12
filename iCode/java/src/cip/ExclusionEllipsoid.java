/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package cip;

/**
 * Excludes locations based on grid aligend ellisoids.
 * @author Thomas Cullison, Colorado School of Mines.
 * @version 2010.10.29
 */
public class ExclusionEllipsoid implements CIPDistanceMap3
{

  /**
   * Constructs an exclusion ellise for specified tensors.
   * @param r1 exlusion radius in the 1st dimension.
   * @param r2 exlusion radius in the 2st dimension.
   * @param r3 exlusion radius in the 3st dimension.
   */
  public ExclusionEllipsoid(float r1, float r2, float r3) {
	  _r1 = r1;
		_r2 = r2;
		_r3 = r3;
		_d1 = 0.f;
		_d2 = 0.f;
		_d3 = 0.f;
		_pii2 = -1;
		_pii3 = -1;
  }

  /**
   * Sets the location (i1,i2) of this brush.
   * @param i1 index in 1st dimension.
   * @param i2 index in 2nd dimension.
   * @param i3 index in 3nd dimension.
   */
  public void setLocation(int i1, int i2, int i3) {
	  _i1=i1;
	  _i2=i2;
	  _i3=i3;
  }


  /**
   * Gets the distance d that the point (ii1,ii2,ii3) is 
	 * from the current location.
   * @param ii1 index in 1st dimension.
   * @param ii2 index in 2nd dimension.
   * @param ii3 index in 3nd dimension.
   */
  public float getDistance(int ii1, int ii2, int ii3) {
		if(ii3 != _pii3) {
			_pii3 = ii3;
			_d3 = _i3-ii3;
			_d3 *= _d3;
			_d3 /= _r3*_r3;
		}
		if(ii2 != _pii2) {
			_pii2 = ii2;
			_d2 = _i2-ii2;
			_d2 *= _d2;
			_d2 /= _r2*_r2;
		}
		_d1 = _i1-ii1;
		_d1 *= _d1;
		_d1 /= _r1*_r1;
		return _d3 + _d2 + _d1;
  }

  private float _r1,_r2,_r3,_d3,_d2,_d1;
  private int _i1,_i2,_i3;
	private int _pii2,_pii3;
}
