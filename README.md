Temporal Plugins
================

This Maven project implements ImageJ 1.x plugins for time sequences:

* Temporal median: finds moving foreground features, see 
    Parton et al. (2011), JCB 194 (1): 121.
* Trails: does simple averaging over a time window, making tracks visible.

Free software, released under the GNU General Public License,
http://www.gnu.org/licenses/gpl.html

Copyright Graeme Ball (2013), graemeball@googlemail.com,
written while working at Micron Oxford: www.micron.ox.ac.uk

The latest .jar files can be found on the [Micron Oxford Website](http://www.micron.ox.ac.uk/microngroup/software/Temporal_plugins.jar)

The maven project structure is derived from:
  https://github.com/imagej/minimal-ij1-plugin
  
Temporal Median filter
----------------------

* set time window, and standard deviations above background for foreground
* time window should be more than 2x larger than time taken for a feature
    to traverse a pixel (NB. total window is 2x half-width +1)
* moving foreground identified by intensity increase relative to background
    average (i.e. median) for a pixel over a given time window
* "soft" segmenation, yielding foreground probability related to excess 
    intensity (in standard deviations) over background level
* crude Anscombe transform applied to data to stabilize the variance
