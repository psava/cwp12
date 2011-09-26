//  ImageLSFFactory.java
//
//  THOMAS CULLISON
//  Colorado School of Mines
//  Center for Wave Phenomena

package cip;

import java.util.*;
import static java.lang.Math.sqrt;
import java.io.IOException;
import edu.mines.jtk.dsp.*;

public class ImageLSFFactory
{

	public static LocalSemblanceFilter getLocalSemblanceFilter(int av, int au) {
	 return new LocalSemblanceFilter(av,au);
  }
  
  public static float[][] getLinearSemblance2D(LocalOrientFilter lof,
																							 LocalSemblanceFilter lsf,
																							 float[][] img) {
		EigenTensors2 d = lof.applyForTensors(img);
    return lsf.semblance(LocalSemblanceFilter.Direction2.V, d, img);
  }
  
  public static float[][][] getLinearSemblance3D(LocalOrientFilter lof,
																								 LocalSemblanceFilter lsf,
																								 float[][][] img) {
		EigenTensors3 d = lof.applyForTensors(img);
    return lsf.semblance(LocalSemblanceFilter.Direction3.W, d, img);
  }
  
  public static float[][][] getPlanarSemblance3D(LocalOrientFilter lof,
																								 LocalSemblanceFilter lsf,
																								 float[][][] img) {
		EigenTensors3 d = lof.applyForTensors(img);
    return lsf.semblance(LocalSemblanceFilter.Direction3.VW, d, img);
  }

}
