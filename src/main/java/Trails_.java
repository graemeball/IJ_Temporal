/*
 * Copyright (c) 2013, Graeme Ball
 * Micron Oxford, University of Oxford, Department of Biochemistry.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/ .
 */

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

/**
 * Trail/average intensities over a time window for an image sequence.
 *
 * @author graemeball@googlemail.com
 */
public class Trails_ implements PlugInFilter {
	protected ImagePlus image;

	// image properties
	private int width;
	private int height;
	private int nt;
	private int nz;
	private int nc;

	// plugin parameters
	public int twh;     // time window half-width for trails

	/**
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
	 */
	@Override
	public int setup(String arg, ImagePlus imp) {
		if (arg.equals("about")) {
			showAbout();
			return DONE;
		}
		image = imp;
		if (image.isHyperStack()) {
    		nt = imp.getNFrames();
    		nz = imp.getNSlices();
    		nc = imp.getNChannels();
		} else {
		    // assume simple stack is a time sequence
		    nt = imp.getStackSize();
		    nz = 1;
		    nc = 1;
		    imp.setDimensions(nc, nz, nt);
		}
		return DOES_8G | DOES_16 | DOES_32 | DOES_RGB
		        | CONVERT_TO_FLOAT | STACK_REQUIRED | NO_CHANGES;
	}

	/**
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	@Override
	public void run(ImageProcessor ip) {
		// get width and height
		width = ip.getWidth();
		height = ip.getHeight();

		if (showDialog()) {
		    if (nt > (2*twh + 1)) {
    			ImagePlus imResult = process(image);
    			imResult.setDimensions(nc, nz, nt);
                imResult.setOpenAsHyperStack(true);
    			imResult.updateAndDraw();
    			imResult.show();
		    } else {
                IJ.showMessage("Insufficient time points, " + nt);
            }
		}
	}

	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Trails");

		gd.addNumericField("time window half-width", 2, 0);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		twh = (int)gd.getNextNumber();

		return true;
	}

	/**
	 * Process each image slice, returning new foreground ImagePlus.
	 * Builds array of pixel arrays for sliding window of time frames.
	 *
	 * @param image (multi-dimensional, i.e. multiple frames)
	 */
	public ImagePlus process(ImagePlus image) {
	    ImageStack inStack = image.getStack();
	    ImageStack stackResult = new ImageStack(width, height);
	    
	    // for all channels and slices, process sliding time window
	    for (int c = 1; c <= nc; c++) {
	        for (int z = 1; z <= nz; z++) {
	            // build initial time window array of pixel arrays
	            float[][] tWinPix = new float[2*twh + 1][width*height];
	            int wmin = 0;  // window min index
	            int wcurr = 0;  // index within window of current frame
	            int wmax = twh;  // window max index
	            for (int t = 1; t <= wmax + 1; t++) {
	                int index = image.getStackIndex(c, z, t);
	                tWinPix[t-1] = getfPixels(inStack, index);
	            }
	            // process each t and update sliding time window
	            for (int t = 1; t <= nt; t++) {
	                //IJ.log("t,wmin,wcurr,wmax=" + t + ","
	                //        + wmin + "," + wcurr + "," + wmax);
	                float[] fgPix = trail(tWinPix, wcurr, wmin, wmax);
	                FloatProcessor fp2 = 
	                        new FloatProcessor(width, height, fgPix);
	                stackResult.addSlice((ImageProcessor)fp2);
	                // sliding window update for next t
	                if (t > twh) {
	                    // remove old pixel array from start
	                    tWinPix = rmFirst(tWinPix, wmax);
	                } else {
	                    wcurr += 1;
	                    wmax += 1;	                    
	                }
	                if (t < nt - twh) {
	                    // append new pixel array (frame t+twh) to end
	                    int newPixIndex = image.getStackIndex(c, z, t+twh);
	                    tWinPix[wmax] = getfPixels(inStack, newPixIndex);
	                } else {
	                    wmax -= 1;	                    
	                }
	            }
	        }
	    }
		return new ImagePlus("Trail" + Integer.toString(2*twh +1) + 
		        "_" + image.getTitle(), stackResult);
	}
	
	/** Remove first array of pixels and shift the others to the left. */
	private float[][] rmFirst(float[][] tWinPix, int wmax) {
	    for (int i=0; i < wmax; i++) {
	        tWinPix[i] = tWinPix[i+1];
	    }
	    return tWinPix;
	}
	
	/** Trail tCurr pixels using tWinPix time window. */
	private float[] trail(float[][] tWinPix, int wcurr, int wmin, int wmax) {
	    int numPix = width*height;
	    float[] tPix = new float[numPix];
	    for (int v=0; v<numPix; v++) {
	        float[] tvec = getTvec(tWinPix, v, wmin, wmax);
	        tPix[v] = fmean(tvec);
	    }
	    return tPix;
	}
	
	/** Build time vector for this pixel for  given window. */
	private float[] getTvec(float[][] tWinPix, int v, int wmin, int wmax) {
	    float[] tvec = new float[wmax - wmin + 1];
	    for (int w=wmin; w<=wmax; w++) {
	        tvec[w] = tWinPix[w][v];  // time window vector for a pixel
	    }
	    return tvec;
	}
	
	/** Calculate mean of array of floats. Shocking. */
	private float fmean(float[] tvec) {
	    float mean = 0;
	    for (int t=0; t<tvec.length; t++) {
	        mean += tvec[t];
	    }
	    return mean / tvec.length;
	}
	
	/** 
	 * Return a float array of pixels for a given stack slice. 
	 */
	private float[] getfPixels(ImageStack stack, int index) {
	    ImageProcessor ip = stack.getProcessor(index);
	    FloatProcessor fp = (FloatProcessor)ip.convertToFloat();
	    float[] pix = (float[])fp.getPixels();
	    return pix;
	}
	
	public void showAbout() {
		IJ.showMessage("Trails",
			"Trail/average intensities over a given time window."
		);
	}

	/**
	 * Main method for debugging. FIXME.
	 * Same as TemporalMedian, but does not work. Sigh.
	 *
	 * For debugging - start ImageJ, load a test image, call the plugin.
	 *
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Trails_.class;
		//String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		//String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		//System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open TrackMate FakeTracks test data from the Fiji wiki
		ImagePlus image = IJ.openImage("http://fiji.sc/tinevez/TrackMate/FakeTracks.tif");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
