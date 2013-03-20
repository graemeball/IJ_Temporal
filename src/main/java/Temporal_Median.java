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
 * A "probabilistic" temporal median filter to extract a foreground
 * probability map from a time sequence.
 *
 * @author graemeball@googlemail.com
 */
public class Temporal_Median implements PlugInFilter {
	protected ImagePlus image;

	// image properties
	private int width;
	private int height;
	private int nt;
	private int nz;
	private int nc;

	// plugin parameters
	public int twh;     // time window half-width for median calc
	public double nsd;  // number of stdev's above median for foreground
	public int value = 1;   // TODO, bin this junk

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
		nt = imp.getNFrames();
		nz = imp.getNSlices();
		nc = imp.getNChannels();
		//return DOES_8G | DOES_16 | DOES_32 | DOES_RGB; 
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
			ImagePlus imResult = process(image);
			imResult.updateAndDraw();
			imResult.show();
		}
	}

	private boolean showDialog() {
		GenericDialog gd = new GenericDialog("Temporal Median");

		// default value is 0.00, 2 digits right of the decimal point
		gd.addNumericField("time window half-width", 5, 0);
		gd.addNumericField("stdevs over median for foreground", 3, 0);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		// get entered values
		twh = (int)gd.getNextNumber();
		nsd = gd.getNextNumber();

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
	                float[] fgPix = calcFg(tWinPix, wcurr, wmin, wmax);
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
		return new ImagePlus("TMFilt_" + image.getTitle(), stackResult);
	}
	
	/** Remove first array of pixels and shift the others to the left. */
	private float[][] rmFirst(float[][] tWinPix, int wmax) {
	    for (int i=0; i < wmax; i++) {
	        tWinPix[i] = tWinPix[i+1];
	    }
	    return tWinPix;
	}
	
	/** Calculate foreground pixel array using for tCurr using tWinPix. */
	private float[] calcFg(float[][] tWinPix, int wcurr, int wmin, int wmax) {
	    return tWinPix[wcurr];
	}
	
	
	/** Return a float array of pixels for a given stack slice. */
	private float[] getfPixels(ImageStack stack, int index) {
	    ImageProcessor ip = stack.getProcessor(index);
	    FloatProcessor fp = (FloatProcessor)ip.convertToFloat();
	    float[] pix = (float[])fp.getPixels();
	    return pix;
	}
	
	public void showAbout() {
		IJ.showMessage("TemporalMedian",
			"A probabilistic temporal median filter, as described in " +
			"Parton et al. (2011), JCB 194 (1): 121."
		);
	}

	/**
	 * Main method for debugging. FIXME.
	 *
	 * For debugging - start ImageJ, load a test image, call the plugin.
	 *
	 * @param args unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins menu
		Class<?> clazz = Temporal_Median.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open TrackMate FakeTracks test data from the Fiji wiki
		ImagePlus image = IJ.openImage("http://fiji.sc/tinevez/TrackMate/FakeTracks.tif");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
