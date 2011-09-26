//  ImageLOFFactory.java
//
//  THOMAS CULLISON
//  Colorado School of Mines
//  Center for Wave Phenomena

package cip;

import java.util.*;
import static java.lang.Math.sqrt;
import java.io.IOException;
import edu.mines.jtk.dsp.*;

public class ImageLOFFactory
{

	public static LocalOrientFilter getLocalOrientFilter(double sig1) {
    return new LocalOrientFilter(sig1);
  }

	public static LocalOrientFilter getLocalOrientFilter(double sig1, 
																											 double sig2) {
    return new LocalOrientFilter(sig1,sig2);
  }

	public static LocalOrientFilter getLocalOrientFilter(double sig1, 
																											 double sig2, 
																											 double sig3) {
    return new LocalOrientFilter(sig1,sig2,sig3);
  }

	public static float[][] getLinearity2D(LocalOrientFilter lof, float[][] img) { 
		int n1 = img[0].length;
		int n2 = img.length;
		float[][] lin = new float[n2][n1];
    lof.apply(img,null,null,null,null,null,null,null,lin);
    return lin;
  }

	public static float[][][] getLinearity3D(LocalOrientFilter lof, float[][][] img) { 
		int n1 = img[0][0].length;
		int n2 = img[0].length;
		int n3 = img.length;
		float[][][] lin = new float[n3][n2][n1];
    lof.apply(img,
							null,null,
							null,null,null,null,null,null,
							null,null,null,null,null,null,
							null,lin);
    return lin;
  }

	public static float[][][] getPlanarity3D(LocalOrientFilter lof, float[][][] img) { 
		int n1 = img[0][0].length;
		int n2 = img[0].length;
		int n3 = img.length;
		float[][][] plane = new float[n3][n2][n1];
    lof.apply(img,
							null,null,
							null,null,null,null,null,null,
							null,null,null,null,null,null,
							plane,null);
    return plane;
  }

	public static float[][][] getNormalAll2D(LocalOrientFilter lof, float[][] img) { 
		int n1 = img[0].length;
		int n2 = img.length;
		float[][][] norm = new float[2][n2][n1];
		lof.applyForNormal(img,norm[0],norm[1]);
		return norm;
  }

	public static float[][][][] getNormalAll3D(LocalOrientFilter lof, float[][][] img) { 
		int n1 = img[0][0].length;
		int n2 = img[0].length;
		int n3 = img.length;
		float[][][][] norm = new float[3][n3][n2][n1];
		lof.applyForNormal(img,norm[0],norm[1],norm[2]);
		return norm;
  }

	public static float[][] getNormalByCoord2D(float[][] coord,
																						 float[][][] normA,
																						 float ox1, float ox2,
																						 float dx1, float dx2) { 
		float s1 = 1.f/dx1;
		float s2 = 1.f/dx2;
		int ni = coord[0].length;
		float[][] normI = new float[2][ni];
		for(int ip=0; ip<ni; ++ip) {
			int i1 = (int) ((coord[0][ip] - ox1)*s1 + 0.5f);
			assert i1 < normA[0][0].length;
			int i2 = (int) ((coord[1][ip] - ox2)*s2 + 0.5f);
			assert i2 < normA[0].length;
			normI[0][ip] = normA[0][i2][i1];
			normI[1][ip] = normA[1][i2][i1];
		}
		return normI;
  }

	public static float[][] getNormalByCoord3D(float[][] coord,
																						 float[][][][] normA,
																						 float ox1, float ox2, float ox3,
																						 float dx1, float dx2, float dx3) { 
		float s1 = 1.f/dx1;
		float s2 = 1.f/dx2;
		float s3 = 1.f/dx3;
		int ni = coord[0].length;
		float[][] normI = new float[3][ni];
		for(int ip=0; ip<ni; ++ip) {
			int i1 = (int) ((coord[0][ip] - ox1)*s1 + 0.5f);
			assert i1 < normA[0][0][0].length;
			int i2 = (int) ((coord[1][ip] - ox2)*s2 + 0.5f);
			assert i2 < normA[0][0].length;
			int i3 = (int) ((coord[2][ip] - ox3)*s3 + 0.5f);
			assert i3 < normA[0].length;
			normI[0][ip] = normA[0][i3][i2][i1];
			normI[1][ip] = normA[1][i3][i2][i1];
			normI[2][ip] = normA[2][i3][i2][i1];
		}
		return normI;
  }
}

