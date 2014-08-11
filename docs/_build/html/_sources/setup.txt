Preliminary Setup Procedures
****************************

Inputs
------

1. Establish a directory that contains data for the specific target for which
   you wish to measure a lens model.  This is the directory from which you will
   run the software.  I call this "uvfit00" for the first run on a given
   dataset, "uvfit01" for the second, etc.

2. Inside this directory, you must ensure the following files are present:

 * config.yaml: This is the configuration file that describes where the
   source of interest is located, what type of model to use for the lens and
   source, the name of the image of the target from your interferometric
   data, the name of the uvfits files containing the interferometric
   visibilities, and a few important processing options as well.  Syntax is
   `yaml <http://www.yaml.org>`_.

 * Image of the target from your interferometric data in fits format.  This
   image is used to set the spatial resolution of the model image (modified by
   an optional oversampling parameter).

 * interferometric visibilities in uvfits format.  

Outputs
-------

posteriorpdf.fits: model parameters for every MCMC iteration, in fits format.
Use astropy (either the fits or table module) to inspect the results directly.

