uvmcmcfit
=========

Detailed documentation is available at http://uvmcmcfit.readthedocs.org/

Python software to fit models to interferometric data, appropriately accounting for gravitational lensing as needed.

--------------------------
 SETUP PROCEDURES

 1. Establish a directory that contains data for the specific target for which
 you wish to measure a lens model.  This is the directory from which you will
 run the software.

 I call this "uvfit00" for the first run on a given dataset, "uvfit01" for
 the second, etc.

 2. Inside this directory, you must ensure the following files are present:

 - "config.yaml": This is the configuration file that describes where the source
 of interest is located, what type of model to use for the lens and source, the
 name of the image of the target from your interferometric data, the name of
 the uvfits files containing the interferometric visibilities, and a few
 important processing options as well.  Syntax is yaml.

 - Image of the target from your interferometric data.  The spatial resolution
 of this image (arcseconds per pixel), modified by an optional oversampling
 parameter, defines the spatial resolution in both the unlensed and lensed
 surface brightness maps.

 - interferometric visibilities for every combination of array configuration,
 sideband, and date observed that you want to model.  

 3. More info about the constraints and priors input files.

 - Lenses: The lenses are assumed to have singular isothermal ellipsoid
 profiles.  

 - Sources: Sources are represented by Gaussian profiles.  Source positions are
 always defined relative to the primary lens, unless there is no lens, in which
 case they are defined relative to the emission centroid defined in
 "config.txt."

--------
 OUTPUTS

"posteriorpdf.fits": model parameters for every MCMC iteration, in fits
format.  Use astropy (either the fits or table module) to inspect the results
directly.

