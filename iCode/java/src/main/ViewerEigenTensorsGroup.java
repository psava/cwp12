/****************************************************************************
Copyright (c) 2010, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

/* This class is a modified version of the EigenTensorsGroup class. */

package rsfjviewer;

import java.awt.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.sgl.*;
import java.util.ArrayList;
import static edu.mines.jtk.ogl.Gl.GL_AMBIENT_AND_DIFFUSE;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A group of ellipsoids that display a list of eigentensors.
 * @author Thomas Cullison and Chris Engelsma, Colorado School of Mines.
 * @version 2010.01.05
 */
public class ViewerEigenTensorsGroup extends Group {


  /**
   * Constructs a new eigentensors group with specified samplings.
   * @param s1 sampling in the first dimension.
   * @param s2 sampling in the second dimension.
   * @param s3 sampling in the third dimension.
   */
  public ViewerEigenTensorsGroup(
    Sampling s1, Sampling s2, Sampling s3)
  {
    _sx = s3;
    _sy = s2;
    _sz = s1;
    _emax = 0.f;
    _etn = new EigenTensorsNode();
    this.addChild(_etn);
    setDefaultStates();
  }

  /**
   * Pulls a tensor from an array of tensors at coordinate locations .
   * @param coord an array of locations. The size is etc[3][nlocations]. 
   * @param etc an array of eigenvectors and eigenvalues at coordinate locations.
	 * The size is etc[9][nlocations]. 
   * @param storing true, if storing; false, otherwise.
	 * etc[0] = u1
	 * etc[1] = u2
	 * etc[2] = u3
	 * etc[3] = w1
	 * etc[4] = w2
	 * etc[5] = w3
	 * etc[6] = e1
	 * etc[7] = e2
	 * etc[8] = e3
   */
  public void pullTensor(float coord[][], float etc[][], boolean storing) {
    double dx = _sx.getDelta();
    double dy = _sy.getDelta();
    double dz = _sz.getDelta();

		int nten = etc[0].length;
		float emaxi;
		_emax = 0;
		for(int i=0; i< nten; ++i) {
			emaxi = max(etc[6][i],etc[7][i],etc[8][i]);
			if (_emax<emaxi)
				_emax = emaxi;
		}
    _etiny = 0.0001f*_emax;
    float scale = 1.0f/sqrt(_emax);

		for(int i=0; i< nten; ++i) {
			float uz = etc[0][i],
						uy = etc[1][i], 
						ux = etc[2][i];
			float wz = etc[3][i],
						wy = etc[4][i], 
						wx = etc[5][i];
			float vz = wy*ux-wx*uy, // v = w cross u
						vy = wx*uz-wz*ux,
						vx = wz*uy-wy*uz;
			float eu = etc[6][i], 
						ev = etc[7][i], 
						ew = etc[8][i];
			if (eu<=_etiny) eu = _etiny;
			if (ev<=_etiny) ev = _etiny;
			if (ew<=_etiny) ew = _etiny;
			float su = scale*sqrt(eu);
			float sv = scale*sqrt(ev);
			float sw = scale*sqrt(ew);
			float xc = coord[i][2], yc = coord[i][1], zc = coord[i][0];
			ux *= su*dx; uy *= su*dy; uz *= su*dz;
			vx *= sv*dx; vy *= sv*dy; vz *= sv*dz;
			wx *= sw*dx; wy *= sw*dy; wz *= sw*dz;
			if (storing) elist.add(new Ellipsoid(xc,yc,zc,ux,uy,uz,vx,vy,vz,wx,wy,wz));
			else eTemp = new Ellipsoid(xc,yc,zc,ux,uy,uz,vx,vy,vz,wx,wy,wz);
			_etn.setLists(elist,eTemp);
		}
  }

  /**
   * Clears the temporary ellipsoid (better method?)
   */
  public void clearTempTensor() {
    eTemp = null;
    _etn.setLists(elist,eTemp);
  }

  /**
   * Clears all stored ellipsoids.
   */
  public void clearAll() {
    eTemp = null;
    elist.removeAll(elist);
    _etn.setLists(elist,eTemp);
  }

  /**
   * Sets the size of the ellipsoids.
   * @param size the size of the ellipsoids (in samples).
   */
  public void setEllipsoidSize(float size) {
    _eSize = size;
    update();
  }

  /**
   * Gets the size of the ellipsoids.
   * @return the ellipsoid size.
   */
  public float getEllipsoidSize() {
    return _eSize;
  }

  /**
   * Clears the stored tensor closest to the given coordinates.
   * @param x the x-coordinate.
   * @param y the y-coordinate.
   * @param z the z-coordinate.
   */
  public void clearClosestTensor(double x, double y, double z) {
    float closestDistance = 1e6f;
    if (elist!=null) {
      float distance = closestDistance;
      Ellipsoid closestEllipsoid = null;
      for (Ellipsoid e:elist) {
        distance = (float)sqrt(pow(x-e.xc,2)+pow(y-e.yc,2)+pow(z-e.zc,2));
        if (distance<closestDistance) {
          closestDistance = distance;
          closestEllipsoid = e;
        }
      }
      if (closestDistance<_closestPossible && closestEllipsoid!=null) {
        elist.remove(closestEllipsoid);
      }
    }
    update();
  }

  /**
   * Returns this eigen tensor node.
   * @return this node.
   */
  public Node getNode() {
    return _etn;
  }

  /**
   * Updates the scene graph.
   */
  public void update() {
    dirtyDraw();
  }

  ////////////////////////////////////////////////////////////////////////////
  // private

  /**
   * An ellipsoid.
   * Ellipsoids consist of twelve numbers for proper definition:
   * center coordinates (x,y,z)
   * vectors representing the three major axes directions (x,y,z)x3
   */
  private class Ellipsoid {
    float xc,yc,zc;
    float ux,uy,uz;
    float vx,vy,vz;
    float wx,wy,wz;
    public Ellipsoid(
      float xc, float yc, float zc, 
      float ux, float uy, float uz,
      float vx, float vy, float vz,
      float wx, float wy, float wz)
    {
      this.xc = xc; this.yc = yc; this.zc = zc;
      this.ux = ux; this.uy = uy; this.uz = uz;
      this.vx = vx; this.vy = vy; this.vz = vz;
      this.wx = wx; this.wy = wy; this.wz = wz;
    }
  }

  private EigenTensors3 _et; // Stored eigentensors
  private Sampling _sx,_sy,_sz; // Eigentensor sampling.
  private float _emax; // Maximum eigenvalue.
  private float _etiny; // Smallest possible eigenvalue.
  private float _eSize = 5.0f; // Maximum ellipsoid radius.
  private float _closestPossible = _eSize; // Closest allowable ellipsoid

  private EigenTensorsNode _etn;
  private ArrayList<Ellipsoid> elist = new ArrayList<Ellipsoid>();
  private Ellipsoid eTemp;

  private class EigenTensorsNode extends Node {

    /**
     * Constructs a new eigen tensor node with null lists.
     */
    public EigenTensorsNode() {
    }

    /**
     * Sets the values of the lists to be drawn.
     * @param elist the list of stuck ellipsoids.
     * @param eTemp the current temporary ellipsoid.
     */
    public void setLists(ArrayList<Ellipsoid> elist, Ellipsoid eTemp) {
      this.elist = elist;
      this.eTemp = eTemp;
      dirtyDraw();
    }

    /**
     * Redraws the scene graph.
     * @param dc the draw context.
     */
    protected void draw(DrawContext dc) {
      if (elist!=null) {
        for (Ellipsoid e:elist)
          _eg.draw(e.xc,e.yc,e.zc,
                   _eSize*e.ux,_eSize*e.uy,_eSize*e.uz,
                   _eSize*e.vx,_eSize*e.vy,_eSize*e.vz,
                   _eSize*e.wx,_eSize*e.wy,_eSize*e.wz);
      }
      if (eTemp!=null)
        _eg.draw(eTemp.xc,eTemp.yc,eTemp.zc,
                 _eSize*eTemp.ux,_eSize*eTemp.uy,_eSize*eTemp.uz,
                 _eSize*eTemp.vx,_eSize*eTemp.vy,_eSize*eTemp.vz,
                 _eSize*eTemp.wx,_eSize*eTemp.wy,_eSize*eTemp.wz);
    }

    private ArrayList<Ellipsoid> elist; // The running list of ellipsoids
    private Ellipsoid eTemp; // A temporary ellipsoid.
    private EllipsoidGlyph _eg = new EllipsoidGlyph();
  }

  /**
   * Sets the state of the group.
   * @param color a color.
   * @return the state set.
   */
  private static StateSet defaultStateSet(Color color) {
    StateSet states = new StateSet();
    ColorState cs = new ColorState();
    cs.setColor(color);
    LightModelState lms = new LightModelState();
    lms.setTwoSide(true);
    MaterialState ms = new MaterialState();
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
    ms.setSpecular(Color.WHITE);
    ms.setShininess(100.0f);
    states.add(cs);
    states.add(lms);
    states.add(ms);
    return states;
  }

  /**
   * Sets the default state of the group.
   */
  private void setDefaultStates() {
    //setStates(defaultStateSet(Color.CYAN));
    setStates(defaultStateSet(Color.RED));
  }

  /**
   * Determines the maximum eigenvalue for proper scaling.
   * @return the maximum eigenvalue.
   */
  private float findMaxEigenvalue() {
    int n1 = _et.getN1();
    int n2 = _et.getN2();
    int n3 = _et.getN3();
    float[] e = new float[3];
    float emax = 0.0f;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          _et.getEigenvalues(i1,i2,i3,e);
          float emaxi = max(e[0],e[1],e[2]);
          if (emax<emaxi)
            emax = emaxi;
        }
      }
    }
    _etiny = 0.0001f*emax;
    return emax;
  }
}
