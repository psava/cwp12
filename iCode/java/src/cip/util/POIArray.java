//
//  POIArray.java
//
//  THOMAS CULLISON
//  Colorado School of Mines
//  Center for Wave Phenomena
package cip.util;

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

import cip.*;

public class POIArray  
{
  public static PointOfInterest[] array1D(int n1) {
	  return new PointOfInterest[n1];
	}

  public static PointOfInterest[][] array2D(int n1, int n2) {
	  return new PointOfInterest[n2][n1];
	}

  public static PointOfInterest[][][] array3D(int n1, int n2, int n3) {
	  return new PointOfInterest[n3][n2][n1];
	}

	private POIArray() {
	}
}
