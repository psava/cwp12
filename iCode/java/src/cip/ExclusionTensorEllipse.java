/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package cip;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.Tensors2;
//import fmm.TimeSolver2;

import static edu.mines.jtk.util.MathPlus.max;
import static edu.mines.jtk.util.MathPlus.min;

/**
 * Excludes locations based on Tensors .
 * @author Thomas Cullison and Chris Englesma, Colorado School of Mines.
 * @version 2010.10.29
 */
public class ExclusionTensorEllipse implements CIPDistanceMap2
{

  /**
   * Constructs an exclusion ellise for specified tensors.
   * @param n1 number of tensor samples in 1st dimension.
   * @param n2 number of tensor samples in 2nd dimension.
   * @param et ellipse tensors.
   */
  public ExclusionTensorEllipse(int n1, int n2, Tensors2 et) {
    _n1 = n1;
    _n2 = n2;
    _et = et;
    setSize(0);
  }

  public void setTensors(Tensors2 et) {
    _et = et;
  }

  /**
   * Gets the location (k1,k2) of this brush.
   * @return array {k1,k2} with location.
   */
  public int[] getLocation() {
    int[] k = new int[2];
    k[0] = _k1;
    k[1] = _k2;
    return k;
  }

  /**
   * Sets the location (k1,k2) of this brush.
   * The default location is (0,0).
   * @param k1 sample index in 1st dimension.
   * @param k2 sample index in 2nd dimension.
   */
  public void setLocation(int k1, int k2) {
    if (_k1!=k1 || _k2!=k2) {
      _k1 = k1;
      _k2 = k2;
      _dirty = true;
    }
  }

  /**
   * Gets the size of this paint brush.
   * @return the size.
   */
  public float getSize() {
    return _size;
  }

  /**
   * Sets the size of this paint brush. For an identity painting tensor
   * (for which time equals distance), the brush size equals its radius.
   * A brush with size zero covers only the one sample where it is located.
   * The default size is ten samples.
   * @param size the size.
   */
  public void setSize(float size) {
    _size = max(0.f,size);

    // Maximum time.
    _tmax = 1.f+size;

    // Half of number of samples in array of brush times; modify as necessary.
    int nh = 1;
    while (nh<=size)
      nh += nh;

    // If number of samples in array of brush times has increased,
    // construct new brush tensors, time solver and marching cubes.
    if (_nh!=nh) {
      int nb = 1+2*nh;
      double db = 1.0;
      double fb = -nh;
      _eet = new ExclusionEllipseTensor2();
      _ts = new TimeSolver2(nb,nb,_eet);
      _nb = nb;
      _nh = nh;
    }

    _dirty = true;
  }


  public float getDistance(int i1, int i2) {
		//System.out.println("size = " + _size);
    int ii1 = i1-_k1+_nh;
    int ii2 = i2-_k2+_nh;

    if (ii1<0 || ii1>=_nb) return INFINITY;
    if (ii2<0 || ii2>=_nb) return INFINITY;

    invokeTimeSolver();

    float time = _ts.getTimes()[ii2][ii1];
		//System.out.println("i2,i1 = " + i2 + "," + i1);
		//System.out.println("k2,k1 = " + _k2 + "," + _k1);
    //float time = _ts.getTimes()[_nh][_nh];
    return time;
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  // Excluder tensors are a subset of the painting tensors.
  private class ExclusionEllipseTensor2 implements Tensors2 {
    public void getTensor(int i1, int i2, float[] a) {
      i1 = max(0,min(_n1-1,i1+_k1-_nh));
      i2 = max(0,min(_n2-1,i2+_k2-_nh));
      _et.getTensor(i1,i2,a);
    }
  }

  private int _n1,_n2;
  private Sampling _s1,_s2;
  private Tensors2 _et;
  private int _k1,_k2;
  private int _nb,_nh;
  private float _size,_tmax;
  private ExclusionEllipseTensor2 _eet;
  private TimeSolver2 _ts;
  private boolean _dirty;
  private static final float INFINITY = Float.MAX_VALUE;

  private void invokeTimeSolver() {
    if (_dirty) {
			//System.out.println("invokedTimeSolver");
      _ts.reset();
      _ts.zeroAt(_nh,_nh);
      _dirty = false;
    }
  }
}
