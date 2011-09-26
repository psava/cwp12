package rsfjviewer;

import edu.mines.jtk.sgl.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import java.util.ArrayList;
import java.io.*;
import java.util.*;
import java.awt.image.*;

import static edu.mines.jtk.ogl.Gl.*;
import java.nio.ByteBuffer;
import edu.mines.jtk.util.*;

import javax.swing.JFileChooser;
import edu.mines.jtk.awt.ColorMap;
import rsf.Input;
import rsf.Output;
import rsf.RSF;

import java.awt.*;
import java.awt.event.ActionEvent;
import javax.swing.*;

import javax.imageio.ImageIO;

import edu.mines.jtk.awt.*;
import edu.mines.jtk.dsp.Sampling;

/****************************************************************************
Copyright (c) 2009, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/

enum ColorList {
    JET(0), GRAY(1), RWB(2), PRISM(3);
    
    private int code;
    private ColorList(int c){
        this.code = c;
    }
    public int getCode(){
        return code;
    }
    
    public static ColorList getMatch(int c){
        for (ColorList color : ColorList.values()){
            if (color.getCode() == c) return color;
        } return GRAY;
    }
}

/**
 * A simple frame for 3D graphics modified for RSF viewing... Borrowed from
 * @author Chris Engelsma and Dave Hale, Colorado School of Mines.
 * @version 2009.07.20
 */
class RSFFrame extends PlotFrame {

  /**
   * Constructs a SimpleFrame with default parameters.
   * Axes orientation defaults to x-right, y-top and z-down.
   */
  public RSFFrame() {
    this(null,null);
  }

  /**
   * Constructs a SimpleFrame with specified axes orientation.
   * @param ao the axes orientation.
   */
  public RSFFrame(AxesOrientation ao) {
    this(null,ao);
  }

  /**
   * Constructs a SimpleFrame with the specified world.
   * @param world the world view.
   */
  public RSFFrame(World world) {
    this(world,null);
  }

  /**
   * Constructs a simple frame with the specified world and orientation.
   * @param world the world view.
   * @param ao the axes orientation.
   */
  public RSFFrame(World world, AxesOrientation ao) {
		super(new PlotPanel());
    if (world==null) world = new World();
    if (ao==null) ao = AxesOrientation.XRIGHT_YOUT_ZDOWN;
    _world = world;
    _view = new OrbitView(_world);
    _view.setAxesOrientation(ao);
    _canvas = new ViewCanvas();
    _canvas.setView(_view);
    _canvas.setBackground(Color.WHITE);
    
    _points = new ArrayList<PointGroup>();
    _lines = new ArrayList<LineGroup>();

		_d = null;
		_tpx = null;
		_tpy = null;
		_ipg = null;

		_etc = null;
		_coord = null;

    ModeManager mm = new ModeManager();
    mm.add(_canvas);
    OrbitViewMode ovm = new OrbitViewMode(mm);
    SelectDragMode sdm = new SelectDragMode(mm);


    JPopupMenu.setDefaultLightWeightPopupEnabled(false);
    ToolTipManager.sharedInstance().setLightWeightPopupEnabled(false);

    JMenu fileMenu = new JMenu("File");
    fileMenu.setMnemonic('F');
    Action exitAction = new AbstractAction("Exit") {
      public void actionPerformed(ActionEvent event) {
        System.exit(0);
      }
    };

    Action cubeAction = new AbstractAction("Add Cube") {
        public void actionPerformed(ActionEvent event){
                    String filename = chooseFile(".");
                    if (filename != null) addRSFCube(filename);
            }
        };
    
    Action lineAction = new AbstractAction("Add Line"){
        public void actionPerformed(ActionEvent event){
                    String filename = chooseFile(".");
                    if (filename != null) addRSFLine(filename);
        }
    };
    
    Action pointAction = new AbstractAction("Add Points"){
        public void actionPerformed(ActionEvent event){
                    String filename = chooseFile(".");
                    if (filename != null) addRSFPoint(filename);
        }
    };
    
    Action loadViewAction = new AbstractAction("Load Viewpoint"){
        public void actionPerformed(ActionEvent event){
                String filename = chooseFile(".");
                if (filename != null) loadView(filename);
        }
    };
    
    Action saveViewAction = new AbstractAction("Save Viewpoint"){
        public void actionPerformed(ActionEvent event){
                JFileChooser chooser = new JFileChooser(new File("."));
                int returnVal = chooser.showSaveDialog(new JFrame());
                String filename = null;
                if(returnVal == JFileChooser.APPROVE_OPTION) {
                        filename = chooser.getSelectedFile().getPath();
                        saveView(filename);
                }
        }
    };
    
    Action saveFrameAction = new AbstractAction("Save to PNG"){
        public void actionPerformed(ActionEvent event){
                JFileChooser chooser = new JFileChooser(new File("."));
                int returnVal = chooser.showSaveDialog(new JFrame());
                String filename = null;
                if(returnVal == JFileChooser.APPROVE_OPTION) {
                        filename = chooser.getSelectedFile().getPath();
                        saveFrametoPNG(filename);
                }
        }
    };


    JMenuItem cubeItem = fileMenu.add(cubeAction);
    cubeItem.setMnemonic('C');
    
    JMenuItem lineItem = fileMenu.add(lineAction);
    lineItem.setMnemonic('L');
    
    JMenuItem pointItem = fileMenu.add(pointAction);
    pointItem.setMnemonic('P');

    JMenuItem saveViewItem = fileMenu.add(saveViewAction);
    saveViewItem.setMnemonic('V');
    
    JMenuItem loadViewItem = fileMenu.add(loadViewAction);
    loadViewItem.setMnemonic('I');
    
    JMenuItem saveFrameItem = fileMenu.add(saveFrameAction);
    saveFrameItem.setMnemonic('S');
    
    JMenuItem exitItem = fileMenu.add(exitAction);
    exitItem.setMnemonic('X');

    
    JMenu colorMenu = new JMenu("Color");

    Action jetAction = new AbstractAction("Jet") {
      public void actionPerformed(ActionEvent event) {
        _color = ColorList.JET;
        setColorMap();
      }
    };
    
    Action prismAction = new AbstractAction("Prism") {
      public void actionPerformed(ActionEvent event) {
        _color = ColorList.PRISM;
        setColorMap();
      }
    };
    
     Action grayAction = new AbstractAction("Gray") {
      public void actionPerformed(ActionEvent event) {
         _color = ColorList.GRAY;
        setColorMap();
      }
    };
    
     Action rwbAction = new AbstractAction("Red-White-Blue") {
      public void actionPerformed(ActionEvent event) {
        _color = ColorList.RWB;
        setColorMap();
      }
    };
    
    colorMenu.add(jetAction);
    colorMenu.add(prismAction);
    colorMenu.add(grayAction);
    colorMenu.add(rwbAction);
    
    
    JMenu clipMenu = new JMenu("% Clip");
    Action clipUp = new AbstractAction("Set max pclip"){
        public void actionPerformed(ActionEvent event){
            String value = JOptionPane.showInputDialog(new JFrame(),
                "Percentile Clip Max (0-100.0):", _pmax);
            try {
                _pmax = Float.parseFloat(value);
                if (_pmax > 100.0f) _pmax = 100.0f;
                if (_pmax < _pmin) _pmax = _pmin+1.0f;
                System.out.printf("pclip: (%f,%f) \n",_pmin,_pmax);
                _ipg.setPercentiles(_pmin,_pmax);   
                  
            } catch (Exception e) {
                System.out.println(e);
            }      
        }
    };
    Action clipDown = new AbstractAction("Set min pclip"){
        public void actionPerformed(ActionEvent event){
             String value = JOptionPane.showInputDialog(new JFrame(),
                String.format("Percentile Clip Min (0-%f):",_pmax), _pmin);
            try {
                _pmin = Float.parseFloat(value);
                if (_pmin < 0.0f) _pmin = 0.0f;
                if (_pmin > _pmax) _pmin = _pmax-1.0f;
                System.out.printf("pclip: (%f,%f) \n",_pmin,_pmax);
                _ipg.setPercentiles(_pmin,_pmax);   
                  
            } catch (Exception e) {
                System.out.println(e);
            }                
        }
    };
    
    clipMenu.add(clipUp);
    clipMenu.add(clipDown);
    

    JMenu modeMenu = new JMenu("Mode");
    modeMenu.setMnemonic('M');
    JMenuItem ovmItem = new JMenuItem(ovm);
    modeMenu.add(ovmItem);
    JMenuItem sdmItem = new JMenuItem(sdm);
    modeMenu.add(sdmItem);

    JMenu tensorMenu = new JMenu("Tensor");

    Action tenLoadAction = new AbstractAction("Load Tensors") {
        public void actionPerformed(ActionEvent event){
                    String filename = chooseFile(".");
                    if (filename != null) loadTensors(filename);
            }
        };

    Action tenCoordLoadAction = new AbstractAction("Load Tensor Coordinates") {
			public void actionPerformed(ActionEvent event){
				String filename = chooseFile(".");
				if (filename != null) loadTensorCoords(filename);
			}
		};

    Action showTenAction = new AbstractAction("Show Tensors at Coordinates") {
			public void actionPerformed(ActionEvent event){

				String filename;
				int sel;

				if(_ipg == null) {
					sel = JOptionPane.showConfirmDialog(null,
										"An RSF Cube (Image) needs to be loaded. "
										+ " Would you like to load one?",
										"Notice",
										JOptionPane.OK_CANCEL_OPTION,
										JOptionPane.QUESTION_MESSAGE);
					if(sel == JOptionPane.OK_OPTION) {
						filename = chooseFile(".");
						if (filename != null) addRSFCube(filename);
					} else {
						JOptionPane.showConfirmDialog(null,
								"Cannot load tensors because an RFC cube was not loaded.",
								"Error",
								JOptionPane.DEFAULT_OPTION,
								JOptionPane.ERROR_MESSAGE);
						return;
					}
				}
				if(_coord == null) {
					sel = JOptionPane.showConfirmDialog(null,
										"The tensors coordinates need to be loaded. "
										+ " Would you like to load them?",
										"Notice",
										JOptionPane.OK_CANCEL_OPTION,
										JOptionPane.QUESTION_MESSAGE);
					if(sel == JOptionPane.OK_OPTION) {
						filename = chooseFile(".");
						if (filename != null) loadCoord(filename);
					} else {
						JOptionPane.showConfirmDialog(null,
								"Error loading coordinates.",
								"Error",
								JOptionPane.DEFAULT_OPTION,
								JOptionPane.ERROR_MESSAGE);
						return;
					}
				} 
				if(_etc == null) {
					sel = JOptionPane.showConfirmDialog(null,
										"The tensors at coordinates need to be loaded. "
										+ " Would you like to load them?",
										"Notice",
										JOptionPane.OK_CANCEL_OPTION,
										JOptionPane.QUESTION_MESSAGE);
					if(sel == JOptionPane.OK_OPTION) {
						filename = chooseFile(".");
						if (filename != null) loadTensorCoords(filename);
						return;
					} else {
						JOptionPane.showConfirmDialog(null,
								"Cannot show tensors because the tensors at coordinates"
								+ " were not loaded.",
								"Error",
								JOptionPane.DEFAULT_OPTION,
								JOptionPane.ERROR_MESSAGE);
						return;
					}
				}
				/*
				if(_d == null) {
					sel = JOptionPane.showConfirmDialog(null,
										"The tensors need to be loaded. "
										+ " Would you like to load them?",
										"Notice",
										JOptionPane.OK_CANCEL_OPTION,
										JOptionPane.QUESTION_MESSAGE);
					if(sel == JOptionPane.OK_OPTION) {
						filename = chooseFile(".");
						if (filename != null) loadTensors(filename);
					} else {
						JOptionPane.showConfirmDialog(null,
								"Cannot show tensors because the tensors were not loaded.",
								"Error",
								JOptionPane.DEFAULT_OPTION,
								JOptionPane.ERROR_MESSAGE);
						return;
					}
				} 
				if(_etg == null) {
					sel = JOptionPane.showConfirmDialog(null,
										"The tensor coordinates need to be loaded. "
										+ " Would you like to load them?",
										"Notice",
										JOptionPane.OK_CANCEL_OPTION,
										JOptionPane.QUESTION_MESSAGE);
					if(sel == JOptionPane.OK_OPTION) {
						filename = chooseFile(".");
						if (filename != null) loadTensorCoords(filename);
						return;
					} else {
						JOptionPane.showConfirmDialog(null,
								"Cannot show tensors because the tensor coordinates"
								+ " were not loaded.",
								"Error",
								JOptionPane.DEFAULT_OPTION,
								JOptionPane.ERROR_MESSAGE);
						return;
					}
				}
				*/
				_world.addChild(_etg);
			}
		};

    Action hideTenAction = new AbstractAction("Hide Tensors at Coordinates") {
        public void actionPerformed(ActionEvent event){
							if(_etg != null)
								_world.removeChild(_etg);
            }
        };

    Action showTenPanAction = new AbstractAction("Show Tensor Panels") {
			public void actionPerformed(ActionEvent event){

				String filename;
				int sel;
				if(_ipg == null) {
					sel = JOptionPane.showConfirmDialog(null,
										"An RSF Cube (Image) needs to be loaded. "
										+ " Would you like to load one?",
										"Notice",
										JOptionPane.OK_CANCEL_OPTION,
										JOptionPane.QUESTION_MESSAGE);
					if(sel == JOptionPane.OK_OPTION) {
						filename = chooseFile(".");
						if (filename != null) addRSFCube(filename);
					} else {
						JOptionPane.showConfirmDialog(null,
								"Cannot load tensors because an RFC cube was not loaded.",
								"Error",
								JOptionPane.DEFAULT_OPTION,
								JOptionPane.ERROR_MESSAGE);
						return;
					}
				}

				if(_d == null) {
					sel = JOptionPane.showConfirmDialog(null,
										"The tensors need to be loaded. "
										+ " Would you like to load them?",
										"Notice",
										JOptionPane.OK_CANCEL_OPTION,
										JOptionPane.QUESTION_MESSAGE);
					if(sel == JOptionPane.OK_OPTION) {
						filename = chooseFile(".");
						if (filename != null) loadTensors(filename);
						addRSFTensorEllipsoids();
					} else {
						JOptionPane.showConfirmDialog(null,
								"Cannot show tensors because the tensors were not loaded.",
								"Error",
								JOptionPane.DEFAULT_OPTION,
								JOptionPane.ERROR_MESSAGE);
						return;
					}
				} else if(_tpx == null || _tpy == null) {
					addRSFTensorEllipsoids();
				} else {
					ImagePanel ipx = _ipg.getImagePanel(Axis.X);
					ImagePanel ipy = _ipg.getImagePanel(Axis.Y);
					ipx.getFrame().addChild(_tpx);
					ipy.getFrame().addChild(_tpy);
				}
			}
		};

    Action hideTenPanAction = new AbstractAction("Hide Tensor Panels") {
        public void actionPerformed(ActionEvent event){
							if(_tpx != null) {
								ImagePanel ipx = _ipg.getImagePanel(Axis.X);
								ipx.getFrame().removeChild(_tpx);
							}
							if(_tpy != null) {
								ImagePanel ipy = _ipg.getImagePanel(Axis.Y);
								ipy.getFrame().removeChild(_tpy);
							}
            }
        };

    JMenuItem tensorItem0 = tensorMenu.add(tenLoadAction);
    tensorItem0.setMnemonic('L');

    JMenuItem tensorItem2 = tensorMenu.add(tenCoordLoadAction);
    tensorItem2.setMnemonic('C');

    JMenuItem tensorItem1 = tensorMenu.add(showTenPanAction);
    tensorItem1.setMnemonic('T');

    JMenuItem tensorItem5 = tensorMenu.add(hideTenPanAction);
    tensorItem5.setMnemonic('R');

    JMenuItem tensorItem4 = tensorMenu.add(showTenAction);
    tensorItem4.setMnemonic('S');

    JMenuItem tensorItem3 = tensorMenu.add(hideTenAction);
    tensorItem3.setMnemonic('H');


    JMenuBar menuBar = new JMenuBar();
    menuBar.add(fileMenu);
    menuBar.add(modeMenu);
    menuBar.add(tensorMenu);
    menuBar.add(colorMenu);
    menuBar.add(clipMenu);

    JToolBar toolBar = new JToolBar(SwingConstants.VERTICAL);
    toolBar.setRollover(true);
    JToggleButton ovmButton = new ModeToggleButton(ovm);
    toolBar.add(ovmButton);
    JToggleButton sdmButton = new ModeToggleButton(sdm);
    toolBar.add(sdmButton);

		_cb = new ColorBar();
		_cb.setWidthMinimum(45);
		_cb.setFont(_cb.getFont().deriveFont(18.f));
		//_ipg.addColorMapListener(_cb);

    ovm.setActive(true);

    this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    this.setSize(new Dimension(SIZE,SIZE));
    this.add(_canvas,BorderLayout.CENTER);
    this.add(toolBar,BorderLayout.WEST);
		this.add(_cb,BorderLayout.EAST);
    this.setJMenuBar(menuBar);
    this.setVisible(true);
  }
  
  
  public void loadView(String filename){
      Point3 point;
      Vector3 tvec;
      double radius;
      double azimuth;
      double elevation;
      double scale;
      double x,y,z;
      double vx,vy,vz;
     
      try {
        if (_ipg == null) throw new Exception("Must load a cube first!");
        Scanner s = new Scanner(new File(filename));
        radius = s.nextDouble();
        x = s.nextDouble();
        y = s.nextDouble();
        z = s.nextDouble();
        point = new Point3(x,y,z);
        azimuth = s.nextDouble();
        elevation = s.nextDouble();
        scale = s.nextDouble();
        vx = s.nextDouble();
        vy = s.nextDouble();
        vz = s.nextDouble();
				tvec = new Vector3(vx,vy,vz);
        Iterator<ImagePanel> itr = _ipg.getImagePanels();
        while (itr.hasNext()){
            ImagePanel ip = itr.next();
            AxisAlignedFrame aaf = ip.getFrame();
            double lx = s.nextDouble();
            double ly = s.nextDouble();
            double lz = s.nextDouble();
            double mx = s.nextDouble();
            double my = s.nextDouble();
            double mz = s.nextDouble();
            
            Point3 min = new Point3(lx,ly,lz);
            Point3 max = new Point3(mx,my,mz);
            aaf.setCorners(min,max);
        }
        _pmax = s.nextFloat();
        int code = s.nextInt();
        _color = ColorList.getMatch(code);
        setColorMap();
        _view.setWorldSphere(new BoundingSphere(point,radius));
				_view.setTranslate(tvec);
        _view.setAzimuth(azimuth);
        _view.setElevation(elevation);
        _view.setScale(scale);
        _ipg.setPercentiles(_pmin,_pmax);
        
      } catch (Exception e){
        System.out.println("Failed to load view point!");
        System.out.println(e);
      }
  }
  
  
  public void saveView(String filename){
    BoundingSphere bs = _ipg.getBoundingSphere(true);
		Vector3 tvec = _view.getTranslate();
    try{
        FileWriter fw = new FileWriter(filename);
        PrintWriter out = new PrintWriter(fw);
        out.println(bs.getRadius());
        Point3 center = bs.getCenter();
        out.printf("%f %f %f \n",center.x,center.y,center.z);
        out.println(_view.getAzimuth());
        out.println(_view.getElevation());
        out.println(_view.getScale());
        out.printf("%f %f %f \n",tvec.x,tvec.y,tvec.z);
        Iterator<ImagePanel> itr = _ipg.getImagePanels();
        while (itr.hasNext()){
            ImagePanel ip = itr.next();
            AxisAlignedFrame aaf = ip.getFrame();
            Point3 min = aaf.getCornerMin();
            Point3 max = aaf.getCornerMax();
            out.printf("%f %f %f %f %f %f\n",min.x,min.y,min.z,max.x,max.y,max.z);           
        }
        out.println(_pmax);
        out.println(_color.getCode());
        out.close();
    } catch (Exception e){
        System.out.println(e);
    }
  
  }
  
  
  public String chooseFile(String path){
    JFileChooser chooser = new JFileChooser(new File(path));
    int returnVal = chooser.showOpenDialog(new JFrame());
    String filename = null;
    if(returnVal == JFileChooser.APPROVE_OPTION) {
            filename = chooser.getSelectedFile().getPath();
    }
    return filename;
  }
  
  public void saveFrametoPNG(String filename){
    final int frameWidth = _canvas.getWidth();
    final int frameHeight = _canvas.getHeight();
    final ByteBuffer pixelsRGB = Direct.newByteBuffer(frameWidth*frameHeight*3);
    _canvas.runWithContext(new Runnable() {
      public void run() {
        //glPushAttrib(GL_PIXEL_MODE_BIT);
        glReadBuffer(GL_BACK);
        glPixelStorei( GL_PACK_ALIGNMENT, 1 );
        glReadPixels( 0, 0, frameWidth, frameHeight, 
	               GL_RGB, GL_UNSIGNED_BYTE, 
	               pixelsRGB );
	    //glPopAttrib();      
	    }
    });     
	int[] pixelInts = new int[ frameWidth * frameHeight ];
    int p = frameWidth * frameHeight * 3; 
    int q; // Index into ByteBuffer
    int i = 0; // Index into target int[]
    int w3 = frameWidth * 3; // Number of bytes in each row
    for (int row = 0; row < frameHeight; row++) {
	    p -= w3;
	    q = p;
	    for (int col = 0; col < frameWidth; col++) {
	      int iR = pixelsRGB.get(q++);
	      int iG = pixelsRGB.get(q++);
	      int iB = pixelsRGB.get(q++);
	      pixelInts[i++] = 
	        0xFF000000 | ((iR & 0x000000FF) << 16) | 
	        ((iG & 0x000000FF) << 8) | (iB & 0x000000FF);
	    }
      }

      // Create a new BufferedImage from the pixeldata.
    BufferedImage bufferedImage = 
	    new BufferedImage( frameWidth, frameHeight, 
			   BufferedImage.TYPE_INT_ARGB);
    bufferedImage.setRGB( 0, 0, frameWidth, frameHeight, 
			    pixelInts, 0, frameWidth );
			    
	try {
	        javax.imageio.ImageIO.write( 
	          bufferedImage, "PNG", new File(filename) );
    } catch (IOException e) {
	        System.out.println( "Error: ImageIO.write." );
	        e.printStackTrace();
    }		    

    /* End code taken from: http://www.felixgers.de/teaching/jogl/imagingProg.html */

    /* 
	  final BufferedImage image = new BufferedImage(
	  this.getWidth(), this.getHeight(), BufferedImage.TYPE_INT_ARGB);
	  Graphics gr = image.getGraphics();
	  this.printAll(gr);
	  gr.dispose();
		try {
			ImageIO.write(image, "PNG", new File(filename));
    } catch (IOException e) {
	        System.out.println( "Error: ImageIO.write." );
	        e.printStackTrace();
    }		    
		*/
  }
  
  public float[] getPointArray(String filename){
        Input file = new Input(filename);
        int n1 = file.getN(1);
        int n2 = file.getN(2);
        assert n1==3;
        
        float[][] data = new float[n2][3];
        file.read(data);
        
        float[] points = new float[n2*3];
        int ii;
        for(int ip = 0; ip < n2; ++ip){
            ii = ip*3;
            points[ii] = data[ip][2];
            points[ii+1] = data[ip][1];
            points[ii+2] = data[ip][0];
        }
        
        return points;
  
  }
  
  public void addRSFLine(String filename){
        System.out.println("add line: "+filename);
        float[] points = getPointArray(filename);
        
        LineGroup lg = new LineGroup(points);
        _lines.add(lg);
        _world.addChild(lg);
  }
  
  public void addRSFPoint(String filename){
        System.out.println("add points: "+filename);
        float[] points = getPointArray(filename);
        
				float[] rgb = new float[points.length];
				for(int i=0; i<points.length; i+=3)
				  rgb[i] = 1.f;
        PointGroup pg = new PointGroup(0.05f,points,rgb);
        _points.add(pg);
        _world.addChild(pg);
  }
  
  public void addRSFCube(String filename){
    System.out.println("add cube: " + filename);
    Input file = new Input(filename);

    int n3 = file.getN(3);
    float d3 = file.getDelta(3);
    float o3 = file.getOrigin(3);
    int n2 = file.getN(2);
    float d2 = file.getDelta(2);
    float o2 = file.getOrigin(2);
    int n1 = file.getN(1);
    float d1 = file.getDelta(1);
    float o1 = file.getOrigin(1);

    Sampling s1 = new Sampling(n1,d1,o1);
    Sampling s2 = new Sampling(n2,d2,o2);
    Sampling s3 = new Sampling(n3,d3,o3);
		_s1 = s1;
		_s2 = s2;
		_s3 = s3;
    float[][][] data = new float[n3][n2][n1];
    file.read(data);
    file.close();
    System.out.println("Done reading data...");
    _ipg = addImagePanels(s1,s2,s3,data);
		_ipg.setClips(ArrayMath.min(data),ArrayMath.max(data));

		_ipg.addColorMapListener(_cb);

    System.out.printf("pclip: (%f,%f) \n", _ipg.getPercentileMin(),_ipg.getPercentileMax());
    _view.setWorldSphere(_ipg.getBoundingSphere(true));
  }

  public void addRSFTensorEllipsoids(){

		_tpx = new TensorsPanel(_s1,_s2,_s3,_d);
		_tpy = new TensorsPanel(_s1,_s2,_s3,_d);
		ImagePanel ipx = _ipg.getImagePanel(Axis.X);
		ImagePanel ipy = _ipg.getImagePanel(Axis.Y);
		ipx.getFrame().addChild(_tpx);
		ipy.getFrame().addChild(_tpy);

		_tpx.setEllipsoidSize(8);
		_tpy.setEllipsoidSize(8);
  }

	public void loadTensors(String filename) {
    System.out.println("add tensor ellipsoids: " + filename);
    Input file = new Input(filename);

		int n1 = _s1.getCount();
		int n2 = _s2.getCount();
		int n3 = _s3.getCount();

    float [][][][] et = new float[9][n3][n2][n1];

    file.read(et);
    file.close();
    System.out.println("Done reading tensors...");

		_d = new EigenTensors3(n1,n2,n3,true);

		float u1,u2,u3;
		float w1,w2,w3;
		for(int i3=0; i3<n3; ++i3)
			for(int i2=0; i2<n2; ++i2)
				for(int i1=0; i1<n1; ++i1) {
					 u1 = et[0][i3][i2][i1];
					 u2 = et[1][i3][i2][i1];
					 u3 = et[2][i3][i2][i1];
					 w1 = et[3][i3][i2][i1];
					 w2 = et[4][i3][i2][i1];
					 w3 = et[5][i3][i2][i1];
					 if(u3 < 0.f) {
					   u1 *= -1.f;
					   u2 *= -1.f;
					   u3 *= -1.f;
					 }
					 _d.setEigenvectorU(i1,i2,i3,u1,u2,u3);
					 if(w3 < 0.f) {
					   w1 *= -1.f;
					   w2 *= -1.f;
					   w3 *= -1.f;
					 }
					 _d.setEigenvectorW(i1,i2,i3,w1,w2,w3);
				}

		_d.setEigenvalues(et[6],et[7],et[8]);
	}

	/*
	public void loadTensorCoords(String filename) {
		System.out.println("add tensor coordinates: "+filename);
		float[] points = getPointArray(filename);
		
		_etg = new EigenTensorsGroup(_s1,_s2,_s3,_d);

		for(int i=0; i<points.length; i+=3)
		  _etg.pullTensor(points[i+0],points[i+1],points[i+2],true);

		_etg.setEllipsoidSize(4.f);
		_world.addChild(_etg);

	}
	*/

	public void loadCoord(String filename) {
		System.out.println("add tensor coordinates: "+filename);
		Input file = new Input(filename);
		int n1 = file.getN(1);
		int n2 = file.getN(2);
		assert n1==3;
		
		_coord = new float[n2][n1];
		file.read(_coord);
	}

	public void loadTensorCoords(String filename) {
		System.out.println("add tensors at coordinates: "+filename);
		Input file = new Input(filename);
		int n1 = file.getN(1);
		int n2 = file.getN(2);
		assert n2==9;
		
		_etc = new float[n2][n1];
		file.read(_etc);

		_etg = new ViewerEigenTensorsGroup(_s1,_s2,_s3);

		_etg.pullTensor(_coord,_etc,true);

		_etg.setEllipsoidSize(4.f);
		_world.addChild(_etg);

	}


  /**
   * Returns a new simple frame with a triangle group.
   * Triangles will be constructed as vertex normals.
   * @param xyz array of packed vertex coordinates.
   * @return a simple frame.
   */
  public static RSFFrame asTriangles(float[] xyz) {
    return asTriangles(true,xyz);
  }

  /**
   * Returns a new simple frame with a triangle group.
   * @param vn true, for vertex normals; false, for triangle normals.
   * @param xyz array of packed vertex coordinates.
   * @return the simple frame.
   */
  public static RSFFrame asTriangles(boolean vn, float[] xyz) {
    return asTriangles(vn,xyz,null);
  }

  /**
   * Returns a new simple frame with a triangle group.
   * Triangles will be constructed as vertex normals.
   * @param xyz array of packed vertex coordinates.
   * @param rgb array of packed color coordinates.
   * @return the simple frame.
   */
  public static RSFFrame asTriangles(float[] xyz, float[] rgb) {
    return asTriangles(true,xyz,rgb);
  }
  /**
   * Returns a new simple frame with a triangle group.
   * @param vn true, for vertex normals; false, for triangle normals.
   * @param sx sampling of x coordinates; may be non-uniform.
   * @param sy sampling of y coordinates; may be non-uniform.
   * @param z array[nx][ny] of z coordinates z = f(x,y).
   */
  public RSFFrame asTriangles(
    boolean vn, Sampling sx, Sampling sy, float[][] z)
  {
    return asTriangles(new TriangleGroup(vn,sx,sy,z));
  }
  
  public void setPercentiles(float min, float max){
    _ipg.setPercentiles(min,max);
  }
  
  public void setColorMap(){
    java.awt.image.IndexColorModel cm = null;
    switch(_color){
        case JET:
            cm = ColorMap.JET;
            break;
        
        case GRAY:
            cm = ColorMap.GRAY;
            break;
        
        case PRISM: 
            cm = ColorMap.PRISM;
            break;
        
        case RWB:
            cm = ColorMap.RED_WHITE_BLUE;
            break;
        
        default: cm = ColorMap.GRAY;
    } 
        _ipg.setColorModel(cm);
  }

  /**
   * Returns a new simple frame with a triangle group.
   * @param vn true, for vertex normals; false, for triangle normals.
   * @param sx sampling of x coordinates; may be non-uniform.
   * @param sy sampling of y coordinates; may be non-uniform.
   * @param z array[nx][ny] of z coordinates z = f(x,y).
   * @param r array[nx][ny] of red color components.
   * @param g array[nx][ny] of green color components.
   * @param b array[nx][ny] of blue color components.
   */
  public RSFFrame asTriangles(
    boolean vn, Sampling sx, Sampling sy, float[][] z,
    float[][] r, float[][] g, float[][] b)
  {
    return asTriangles(new TriangleGroup(vn,sx,sy,z,r,g,b));
  }

  /**
   * Returns a new simple frame with a triangle group.
   * @param vn true, for vertex normals; false, for triangle normals
   * @param xyz array of packed vertex coordinates.
   * @param rgb array of packed color coordinates.
   * @return the simple frame.
   */
  public static RSFFrame asTriangles(boolean vn, float[] xyz, float[] rgb) {
    return asTriangles(new TriangleGroup(vn,xyz,rgb));
  }

  /**
   * Returns a new simple frame with a triangle group.
   * @param tg a triangle group.
   * @return the simple frame.
   */
  public static RSFFrame asTriangles(TriangleGroup tg) {
    RSFFrame sf = new RSFFrame();
    sf.addTriangles(tg);
    sf.getOrbitView().setWorldSphere(tg.getBoundingSphere(true));
    return sf;
  }

  /**
   * Returns a new simple frame with an image panel group.
   * @param f a 3D array.
   * @return the simple frame.
   */
  public static RSFFrame asImagePanels(float[][][] f) {
    return asImagePanels(new ImagePanelGroup(f));
  }

  /**
   * Returns a new simple frame with an image panel group.
   * @param s1 sampling in the 1st dimension (Z).
   * @param s2 sampling in the 2nd dimension (Y).
   * @param s3 sampling in the 3rd dimension (X).
   * @param f a 3D array.
   * @return the simple frame.
   */
  public static RSFFrame asImagePanels(
          Sampling s1, Sampling s2, Sampling s3, float[][][] f) {
    return asImagePanels(new ImagePanelGroup(s1,s2,s3,f));
  }

  /**
   * Returns a new simple frame with an image panel group.
   * @param ipg an image panel group.
   * @return the simple frame.
   */
  public static RSFFrame asImagePanels(ImagePanelGroup ipg) {
    RSFFrame sf = new RSFFrame();
    sf.addImagePanels(ipg);
    sf.getOrbitView().setWorldSphere(ipg.getBoundingSphere(true));
    return sf;
  }

  /**
   * Adds a triangle group with specified vertex coordinates.
   * @param xyz array of packed vertex coordinates.
   * @return the triangle group.
   */
  public TriangleGroup addTriangles(float[] xyz) {
    return addTriangles(new TriangleGroup(true,xyz,null));
  }

  /**
   * Adds a triangle group with specified vertex coordinates and colors.
   * @param xyz array of packed vertex coordinates.
   * @param rgb array of packed color components.
   * @return the triangle group.
   */
  public TriangleGroup addTriangles(float[] xyz, float[] rgb) {
    return addTriangles(new TriangleGroup(true,xyz,rgb));
  }

  /**
   * Adds a triangle group for a sampled function z = f(x,y).
   * @param sx sampling of x coordinates; may be non-uniform.
   * @param sy sampling of y coordinates; may be non-uniform.
   * @param z array[nx][ny] of z coordinates z = f(x,y).
   */
  public TriangleGroup addTriangles(Sampling sx, Sampling sy, float[][] z) {
    return addTriangles(new TriangleGroup(true,sx,sy,z));
  }

  /**
   * Adds a triangle group to a simple frame from a given triangle group.
   * @param tg a triangle group.
   * @return the attached triangle group.
   */
  public TriangleGroup addTriangles(TriangleGroup tg) {
    _world.addChild(tg);
    return tg;
  }

  /**
   * Adds an image panel group to a simple frame from a given 3D array
   * @param f a 3D array.
   * @return the image panel group.
   */
  public ImagePanelGroup addImagePanels(float[][][] f) {
		Sampling ts1 = new Sampling(f[0][0].length);
		Sampling ts2 = new Sampling(f[0].length);
		Sampling ts3 = new Sampling(f.length);
		_s1 = ts1;
		_s2 = ts2;
		_s3 = ts3;
		return addImagePanels(ts1,ts2,ts3,f);
		/*
    return addImagePanels(new Sampling(f[0][0].length),
                          new Sampling(f[0].length),
                          new Sampling(f.length),
                          f);
		*/
  }

  /**
   * Adds an image panel group to a simple frame from given samplings and
   * a 3D array.
   * @param s1 sampling in the 1st dimension (Z).
   * @param s2 sampling in the 2nd dimension (Y).
   * @param s3 sampling in the 3rd dimension (X).
   * @param f a 3D array.
   * @return the image panel group.
   */
  public ImagePanelGroup addImagePanels(
          Sampling s1, Sampling s2, Sampling s3, float[][][] f) {
    return addImagePanels(new ImagePanelGroup(s1,s2,s3,f));
  }

  /**
   * Adds an image panel group to a simple frame from a given image panel
   * group.
   * @param ipg an image panel group.
   * @return the attached image panel group.
   */
  public ImagePanelGroup addImagePanels(ImagePanelGroup ipg) {
    _world.addChild(ipg);
    return ipg;
  }

  /**
   * Gets the view canvas for this frame.
   * @return the view canvas.
   */
  public ViewCanvas getViewCanvas() {
    return _canvas;
  }

  /**
   * Gets the orbit view for this frame.
   * @return the orbit view.
   */
  public OrbitView getOrbitView() {
    return _view;
  }

  /**
   * Gets the world for this frame.
   * @return the world.
   */
  public World getWorld() {
    return _world;
  }

  /**
   * Sets the bounding sphere of the frame with a given center point and
   * radius.
   * @param p the center point.
   * @param r the radius.
   */
  public void setWorldSphere(Point3 p, int r) {
    setWorldSphere(new BoundingSphere(p,r));
  }

  /**
   * Sets the bounding sphere of the frame with a given center x, y, z and
   * radius.
   * @param x the center X-coordinate.
   * @param y the center Y-coordinate.
   * @param z the center Z-coordinate.
   * @param r the radius.
   */
  public void setWorldSphere(double x, double y, double z, double r) {
    setWorldSphere(new BoundingSphere(x,y,z,r));
  }

  /**
   * Sets the bounding sphere of the frame.
   * @param xmin the minimum x coordinate.
   * @param ymin the minimum y coordinate.
   * @param zmin the minimum z coordinate.
   * @param xmax the maximum x coordinate.
   * @param ymax the maximum y coordinate.
   * @param zmax the maximum z coordinate.
   */
  public void setWorldSphere(
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax)
  {
    setWorldSphere(
      new BoundingSphere(
        new BoundingBox(xmin,ymin,zmin,xmax,ymax,zmax)));
  }

  /**
   * Sets the bounding sphere of the frame.
   * @param bs the bounding sphere.
   */
  public void setWorldSphere(BoundingSphere bs) {
    _view.setWorldSphere(bs);
  }

/////////////////////////////////////////////////////////////////////////////
// private

  private ArrayList<PointGroup> _points;
  private ArrayList<LineGroup> _lines;
  private ViewCanvas _canvas;
  private OrbitView _view;
  private World _world;
  private ImagePanelGroup _ipg;
  private static final int SIZE = 600;
  private float _pmin=0.0f, _pmax=100.0f;
  private ColorList _color=ColorList.GRAY;
	private ColorBar _cb;
	private float[][] _coord;
	private float[][] _etc;
	private EigenTensors3 _d;
	private ViewerEigenTensorsGroup _etg;
	private TensorsPanel _tpx,_tpy;
	private Sampling _s1, _s2, _s3;
}


public class Viewer3D {

	static {
		System.loadLibrary("jrsf");
	}

	public static void main(String[] args){
//      args = new String[]{"nothing=junk"};
      RSF par = new RSF(args);
			//RSFFrame sf = new RSFFrame();
            //
            //
            //
            //
     
      ArrayList<String> cubes = new ArrayList<String>();
      ArrayList<String> points = new ArrayList<String>();
      ArrayList<String> lines = new ArrayList<String>();
      for(String arg: args){
          if (arg.contains("cube")){
              String cube = par.getString(arg.split("=")[0],"");
              System.err.printf("Found cube: %s\n",cube);
              cubes.add(cube);
          } else if (arg.contains("point")){
              String point = par.getString(arg.split("=")[0],"");
              System.err.printf("Found point: %s\n",point);
              points.add(point);
          } else if (arg.contains("line")){
              String line = par.getString(arg.split("=")[0],"");
              System.err.printf("Found line: %s\n",line);
              lines.add(line);
          }
      }

      final String[] cubeNames = new String[cubes.size()];
      for(int i = 0; i < cubes.size(); ++i){
          cubeNames[i] = cubes.get(i);
      }
      
      final String[] lineNames = new String[lines.size()];
      for(int i = 0; i < lines.size(); ++i){
          lineNames[i] = lines.get(i);
      }
      final String[] pointNames =  new String[points.size()];
      for(int i = 0; i < points.size(); ++i){
          pointNames[i] = points.get(i);
      }

      SwingUtilities.invokeLater(new Runnable() {
          public void run() {
            RSFFrame frame = new RSFFrame();
            for (String name: cubeNames){
                frame.addRSFCube(name);
            }
            for (String name: lineNames){
                frame.addRSFLine(name);
            }
            for (String name: pointNames){
                frame.addRSFPoint(name);
            }
      }
    });
	}
}
