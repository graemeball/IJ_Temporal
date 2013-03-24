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
import java.util.Arrays;

/**
 * A "probabilistic" temporal median filter to extract a foreground
 * probability image from a time sequence.
 *
 * @author graemeball@googlemail.com
 */
public class TemporalMedian_ implements PlugInFilter {

	protected ImagePlus image;

	// image properties
	private int width;
	private int height;
	private int nt;
	private int nz;
	private int nc;

	// plugin parameters
	public int twh;     // time window half-width for median calc
	public float nsd;  // number of stdev's above median for foreground

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
		GenericDialog gd = new GenericDialog("Temporal Median");

		gd.addNumericField("time window half-width", 5, 0);
		gd.addNumericField("stdevs over median for foreground", 2.0, 1);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;

		// get entered values
		twh = (int)gd.getNextNumber();
		nsd = (float)gd.getNextNumber();
		
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
	    int progressCtr = 0;
	    IJ.showStatus("Finding Foreground...");
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
	                    int newPixIndex = image.getStackIndex(c, z, t+twh+1);
	                    tWinPix[wmax] = getfPixels(inStack, newPixIndex);
	                } else {
	                    wmax -= 1;	                    
	                }
	                IJ.showProgress(progressCtr++, nc*nz*nt);
	            }
	        }
	    }
		return new ImagePlus("TMFilt_" + image.getTitle(), stackResult);
	}
	
	/** 
	 * Return a float array of variance-stabilized pixels for a given 
	 * stack slice - applies Anscombe transform. 
	 */
	private float[] getfPixels(ImageStack stack, int index) {
	    ImageProcessor ip = stack.getProcessor(index);
	    FloatProcessor fp = (FloatProcessor)ip.convertToFloat();
	    float[] pix = (float[])fp.getPixels();
	    for (int i=0; i<pix.length; i++) {
	        double raw = (double)pix[i];
	        double transf = 2*Math.sqrt(raw + 3/8);
	        pix[i] = (float)transf;
	    }
	    return pix;
	}
	
	/** Calculate foreground pixel array using for tCurr using tWinPix. */
	private float[] calcFg(float[][] tWinPix, int wcurr, int wmin, int wmax) {
	    float sd = estimStdev(tWinPix, wmin, wmax);
	    int numPix = width*height;
	    float[] fgPix = new float[numPix];
	    for (int v=0; v<numPix; v++) {
	        float[] tvec = getTvec(tWinPix, v, wmin, wmax);
	        float median = fmedian(tvec);
	        float currPix = tWinPix[wcurr][v]; 
	        fgPix[v] = calcFgProb(currPix, median, sd);
	    }
	    return fgPix;
	}
	
	/** Build time vector for this pixel for  given window. */
	private float[] getTvec(float[][] tWinPix, int v, int wmin, int wmax) {
	    float[] tvec = new float[wmax - wmin + 1];
	    for (int w=wmin; w<=wmax; w++) {
	        tvec[w] = tWinPix[w][v];  // time window vector for a pixel
	    }
	    return tvec;
	}
	
	/** Calculate median of an array of floats. */
	private float fmedian(float[] m) {
	    Arrays.sort(m);
	    // as suggested by Nico Huysamen, S.O. q.4191687
	    int middle = m.length/2;
	    if (m.length % 2 == 1) {
	        return m[middle];
	    } else {
	        return (m[middle-1] + m[middle]) / 2.0f;
	    }
	}
	
	/** Remove first array of pixels and shift the others to the left. */
	private float[][] rmFirst(float[][] tWinPix, int wmax) {
	    for (int i=0; i < wmax; i++) {
	        tWinPix[i] = tWinPix[i+1];
	    }
	    return tWinPix;
	}
	
	/** 
	 * Estimate Stdev for this time window using random 0.1% of tvecs. 
	 * Returns the average (mean) stdev of a sample of random tvecs.
	 */
	private float estimStdev(float[][] tWinPix, int wmin, int wmax) {
	    float sd = 0;
	    int pixArrayLen = tWinPix[0].length;
	    int samples = tWinPix.length * pixArrayLen / 1000;
	    for (int n=0; n<samples; n++) {
	        int randPix = (int)Math.floor(Math.random()*pixArrayLen);
	        float[] tvec = getTvec(tWinPix, randPix, wmin, wmax);
	        sd += calcSD(tvec)/samples;
	    }
	    return sd;
	}
	
	/** Standard deviation of a vector of float values. */
	private float calcSD(float[] vec) {
	    float sd = 0;
	    float mean = 0;
	    float variance = 0;
	    for (float v : vec) {
	        mean += v;
	    }
	    mean /= vec.length;
	    for (float v: vec) {
	        variance += (mean - v) * (mean - v);
	    }
	    variance /= vec.length;
	    sd = (float)Math.sqrt(variance);
	    return sd;
	}
	
	/** 
	 * Calculate foreground probability for a pixel using tvec median & stdev. 
	 * foreground probability, P(x,y,z,t) = Q(-v), where:
	 * v = [I(x,y,z,t) - Ibg(x,y,z,t) - k*sigma]/sigma ;  
	 * I and Ibg are Intensity and background intensity (i.e. temporal median); 
	 * sigma is standard deviation of intensity over time ; 
	 * Q is the Q-function (see calcQ).
	 * 
	 */
	float calcFgProb(float currPix, float median, float sd) {
	    float fgProb;
	    fgProb = (currPix - median - nsd*sd)/sd;
	    fgProb = calcQ(-fgProb);
    	return fgProb;
	}
	
	/**
	 * Calculate Q-function, Q(v) = 0.5*(1 - erf(v/sqrt(2))) ; 
     * where erf is the error function ;
	 * see: see http://en.wikipedia.org/wiki/Q-function
	 */
	private static float calcQ(float v) {
	    float root2 = 1.4142135623730950488016887f; // magic ;-)
	    float Q;
	    Q = (0.5f * (1.0f - calcErf(v / root2)));
	    return Q;
	}
	
	/**
	 * Calculate the error function. See: 
	 * http://en.wikipedia.org/wiki/Error_function
	 */
	private static float calcErf(float v) {
	    // room for improvement - this approximation is fast & loose
	    float erf;
	    float a1 = 0.278393f;
	    float a2 = 0.230389f;
	    float a3 = 0.000972f;
	    float a4 = 0.078108f;
	    boolean neg = false;
	    // to use approximation for -ve values of 'v', use: erf(v) = -erf(-v)
	    if (v < 0) {
	        v = -v;
	        neg = true;
	    }
	    erf = 1.0f - (float)(1.0 / Math.pow((double)(1.0 + 
                                        a1 * v + 
                                        a2 * v * v +
                                        a3 * v * v * v +
                                        a4 * v * v * v * v), 4.0));
	    if (neg) {
	        erf = -erf;
	    }
	    return erf;
	}
	
	public void showAbout() {
		IJ.showMessage("TemporalMedian",
			"A probabilistic temporal median filter, as described in " +
			"Parton et al. (2011), JCB 194 (1): 121."
		);
	}

	/**
	 * Main method for debugging.
	 *
	 * For debugging - start ImageJ, load a test image, call the plugin.
	 *
	 * @param args unused
	 */
	public static void main(String[] args) {
		Class<?> clazz = TemporalMedian_.class;
		
		// print calcErf and calcQ results in range -2->2 to check
		for (float i=-2; i<2; i+=0.2) {
		    System.out.println("erf(" + i + ") = " + calcErf(i));
		    System.out.println("Q(" + i + ") = " + calcQ(i));
		}

		// start ImageJ
		new ImageJ();

		// open TrackMate FakeTracks test data from the Fiji wiki
		ImagePlus image = IJ.openImage("http://fiji.sc/tinevez/TrackMate/FakeTracks.tif");
		image.show();

		// run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}
}
