//  ImageTensorFactory.java
//
//  THOMAS CULLISON
//  Colorado School of Mines
//  Center for Wave Phenomena

package cip;

import java.util.*;
import static java.lang.Math.sqrt;
import java.io.IOException;
import edu.mines.jtk.dsp.*;

public class ImageTensorFactory
{

	public static float[][][] getMetricTensors2D(LocalOrientFilter lof, 
																							 LocalSemblanceFilter lsf,
																							 float[][] img) { 
		int n1 = img[0].length;
		int n2 = img.length;
		float[][][] et = new float[4][n2][n1];
		lof.apply(img,null,
		          et[0],et[1],
							null,null,
							null,null,null);
		EigenTensors2 d = lof.applyForTensors(img);
		et[2] = lsf.semblance(LocalSemblanceFilter.Direction2.UV, d, img);
		et[3] = lsf.semblance(LocalSemblanceFilter.Direction2.V, d, img);
		return et;
  }

	public static float[][][][] getMetricTensors3D(LocalOrientFilter lof, 
																								 LocalSemblanceFilter lsf,
																								 float[][][] img) { 
		int n1 = img[0][0].length;
		int n2 = img[0].length;
		int n3 = img.length;
		float[][][][] et = new float[9][n3][n2][n1];
		lof.apply(img,null, null,
		          et[0],et[1],et[2],
							null,null,null,
							et[3],et[4],et[5],
							null,null,null,null,null);
		EigenTensors3 d = lof.applyForTensors(img,true);
    et[6] = lsf.semblance(LocalSemblanceFilter.Direction3.UVW, d, img);
    et[7] = lsf.semblance(LocalSemblanceFilter.Direction3.VW, d, img);
    et[8] = lsf.semblance(LocalSemblanceFilter.Direction3.W, d, img);
		return et;
  }

	public static float[][][] getStructureTensors2D(LocalOrientFilter lof, 
																									float[][] img) { 
		int n1 = img[0].length;
		int n2 = img.length;
		float[][][] et = new float[4][n2][n1];
		lof.apply(img,null,
		          et[0],et[1],
							null,null,
							null,null,null);
		EigenTensors2 d = lof.applyForTensors(img);
		d.getEigenvalues(et[2], et[3]);
		return et;
  }

	public static float[][][][] getStructureTensors3D(LocalOrientFilter lof, 
																										float[][][] img) { 
		int n1 = img[0][0].length;
		int n2 = img[0].length;
		int n3 = img.length;
		float[][][][] et = new float[9][n3][n2][n1];
		lof.apply(img,null, null,
		          et[0],et[1],et[2],
							null,null,null,
							et[3],et[4],et[5],
							null,null,null,null,null);
		EigenTensors3 d = lof.applyForTensors(img,true);
		d.getEigenvalues(et[6], et[7], et[8]);
		return et;
  }

	public static float[][][] getOrientedTensors2D(LocalOrientFilter lof, 
																								 float[][] img) { 
		int n1 = img[0].length;
		int n2 = img.length;
		float[][][] et = new float[4][n2][n1];
		lof.apply(img,null,
		          et[0],et[1],
							null,null,
							null,null,null);
		EigenTensors2 d = lof.applyForTensors(img);
		d.setEigenvalues(0.05f,1.f);
		d.getEigenvalues(et[2], et[3]);
		return et;
  }

	public static float[][][][] getOrientedTensors3D(LocalOrientFilter lof, 
																									 float[][][] img) { 
		int n1 = img[0][0].length;
		int n2 = img[0].length;
		int n3 = img.length;
		float[][][][] et = new float[9][n3][n2][n1];
		lof.apply(img,null, null,
		          et[0],et[1],et[2],
							null,null,null,
							et[3],et[4],et[5],
							null,null,null,null,null);
		EigenTensors3 d = lof.applyForTensors(img,true);
		d.setEigenvalues(0.05f,1.f,1.f);
		d.getEigenvalues(et[6], et[7], et[8]);
		return et;
  }

	public static float[][][] getExclusionTensors2D(LocalOrientFilter lof, 
																								  float[][] img,
																									float ecc) { 
		int n1 = img[0].length;
		int n2 = img.length;
		float[][][] et = new float[4][n2][n1];
		lof.apply(img,null,
		          et[0],et[1],
							null,null,
							null,null,null);
		EigenTensors2 d = lof.applyForTensors(img);
		float eccU,eccV;
		if(ecc >= 1.f) {
		  eccU = 1.f;
		  eccV = 1.f/ecc;
		} else {
		  eccU = ecc;
		  eccV = 1.f;
		}
		d.setEigenvalues(eccU,eccV);
		d.getEigenvalues(et[2],et[3]);
		return et;
  }

	public static float[][][][] getExclusionTensors3D(LocalOrientFilter lof, 
																									  float[][][] img,
																										float ecc) { 
		int n1 = img[0][0].length;
		int n2 = img[0].length;
		int n3 = img.length;
		float[][][][] et = new float[9][n3][n2][n1];
		lof.apply(img,null, null,
		          et[0],et[1],et[2],
							null,null,null,
							et[3],et[4],et[5],
							null,null,null,null,null);
		EigenTensors3 d = lof.applyForTensors(img,true);
		float eccU,eccVW;
		if(ecc >= 1.f) {
		  eccU = 1.f;
		  eccVW = 1.f/ecc;
		} else {
		  eccU = ecc;
		  eccVW = 1.f;
		}
		d.setEigenvalues(eccU,eccVW,eccVW);
		d.getEigenvalues(et[6],et[7],et[8]);
		return et;
  }

	public static float[][] getTensorsByCoord2D(float[][] coord,
																						 float[][][] eti,
																						 float ox1, float ox2,
																						 float dx1, float dx2) { 
		float s1 = 1.f/dx1;
		float s2 = 1.f/dx2;
		int ni = coord[0].length;
		float[][] eto = new float[4][ni];
		for(int ip=0; ip<ni; ++ip) {
			int i1 = (int) ((coord[0][ip] - ox1)*s1 + 0.5f);
			assert i1 < eti[0][0].length;
			int i2 = (int) ((coord[1][ip] - ox2)*s2 + 0.5f);
			assert i2 < eti[0].length;
			eto[0][ip] = eti[0][i2][i1];
			eto[1][ip] = eti[1][i2][i1];
			eto[2][ip] = eti[2][i2][i1];
			eto[3][ip] = eti[3][i2][i1];
		}
		return eto;
  }

	public static float[][] getTensorsByCoord3D(float[][] coord,
																						 float[][][][] eti,
																						 float ox1, float ox2, float ox3,
																						 float dx1, float dx2, float dx3) { 
		float s1 = 1.f/dx1;
		float s2 = 1.f/dx2;
		float s3 = 1.f/dx3;
		int ni = coord[0].length;
		float[][] eto = new float[9][ni];
		for(int ip=0; ip<ni; ++ip) {
			int i1 = (int) ((coord[0][ip] - ox1)*s1 + 0.5f);
			assert i1 < eti[0][0][0].length;
			int i2 = (int) ((coord[1][ip] - ox2)*s2 + 0.5f);
			assert i2 < eti[0][0].length;
			int i3 = (int) ((coord[2][ip] - ox3)*s3 + 0.5f);
			assert i3 < eti[0].length;
			eto[0][ip] = eti[0][i3][i2][i1];
			eto[1][ip] = eti[1][i3][i2][i1];
			eto[2][ip] = eti[2][i3][i2][i1];
			eto[3][ip] = eti[3][i3][i2][i1];
			eto[4][ip] = eti[4][i3][i2][i1];
			eto[5][ip] = eti[5][i3][i2][i1];
			eto[6][ip] = eti[6][i3][i2][i1];
			eto[7][ip] = eti[7][i3][i2][i1];
			eto[8][ip] = eti[8][i3][i2][i1];
		}
		return eto;
  }

}

