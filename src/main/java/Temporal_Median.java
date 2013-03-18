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

	// image property members
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
	 *
	 * @param image (multi-dimensional, i.e. multiple frames)
	 */
	public ImagePlus process(ImagePlus image) {
	    ImageStack stackResult = new ImageStack(width, height);
	    
		// slice numbers start with 1 for historical reasons
		for (int i = 1; i <= image.getStackSize(); i++) {
		    ImageProcessor ip = image.getStack().getProcessor(i);
		    int t = image.getT();
		    int z = image.getZ();
		    int c = image.getC();
			ImageProcessor pip = process(ip, t, z, c);
			stackResult.addSlice(pip);
		}
		
		return new ImagePlus("TMFilt_" + image.getTitle(), stackResult);
	}
	
	/** Temporal median filter / foreground probability calc for a slice */
	public ImageProcessor process(ImageProcessor ip, int t, int z, int c) {
	    
	    return ip;
	}
/*
	// processing of GRAY32 images
	public void process(float[] pixels) {
		for (int y=0; y < height; y++) {
			for (int x=0; x < width; x++) {
				// process each pixel of the line
				// example: add 'number' to each pixel
				pixels[x + y * width] += (float)value;
			}
		}
	}
*/
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
