//
//  PointOfInterest.java
//
//  THOMAS CULLISON
//  Colorado School of Mines
//  Center for Wave Phenomena

package cip;

import java.util.*;
import java.io.IOException;
import java.nio.ByteOrder;
import javax.swing.*;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.io.ArrayInputStream;
import static edu.mines.jtk.util.ArrayMath.*;

public class PointOfInterest implements Comparable<PointOfInterest>
{

	public PointOfInterest() {
		_val = 0.f;
		_dim = 0;
		_index = null;
		_live = false;
		_lock = false;
	}

	public PointOfInterest(float val, int dim, int[] index) {
		_val = val;
		_dim = dim;
		setIndecies(index);
		_live = true;
		_lock = false;
	}

	public int compareTo(PointOfInterest other) {
		if(_val < other.getVal())      return  1;
		if(_val > other.getVal())      return -1;
		return 0;
	}
	
	public void setVal(float val) {
		_val = val;
	}

	public void setIndex(int ix, int i) {
	  _index[ix] = i;
	}

	public void setLive(boolean live) {
		_live = live;
	}

	public void setLock(boolean lock) {
		_lock = lock;
	}

	public float getVal() {
		return _val;
	}

	public float getDim() {
		return _dim;
	}

	public int getIndex(int ix) {
		return _index[ix];
	}

	public boolean isLive() {
		return _live;
	}

	public boolean isLocked() {
		return _lock;
	}

	private void setIndecies(int[] index) {
		_index = new int[_dim];
	  for (int i=0; i<_dim; ++i)
		  _index[i] = index[i];
	}

	private float _val;
	private int   _dim;
	private int[] _index;
	private boolean _live;
	private boolean _lock;
}
