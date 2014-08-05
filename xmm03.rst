XMM03: a single, unlensed galaxy
********************************

This example shows how to run ``uvmcmcfit`` on the simplest of systems: a
single galaxy that is unaffected by lensing of any kind.

.. toctree::
    :maxdepth: 2

    setup

Configuring config.yaml
-----------------------

config.yaml contains the instructions needed by ``uvmcmcfit`` to initiate the
model fitting process.

Required keywords
^^^^^^^^^^^^^^^^^

The name of the target (used for plotting the best-fit model)::

    ObjectName: XMM03

The name of the fits image of the target::
 
    ImageName: XMM03.concat.statwt.cont.mfs.fits

The name of the uvfits visibility data::

    UVData: XMM03.concat.statwt.cont.uvfits

The number of walkers for emcee to use (must be more than double the number of
parameters).  In this case, there are only 6 parameters, so the minimum number
of walkers is 12.  I selected 24 to be on the safe side::
    
    Nwalkers: 24

You must define at least one region.  The first region should be named
*Region0*, the second *Region1*, etc.  Pay attention to the indentation; the
remaining keywords must be indented to indicate they are sub-components of
*Region0*::

    Region0:

    Right Ascension center of the model image (degrees)::

        RACentroid:

Optional keywords
^^^^^^^^^^^^^^^^^

By default, the maximum likelihood estimate is used to measure the goodness of
fit.  Alternatively, you may use the chi-squared value as the goodness of the
fit via::

    LogLike: chi2

By default, parallel processing is not used.  To use parallel processing on a
single machine, set the Nthreads variable to a number greater than 1.  For
example, ::

    Nthreads: 2

If you have access to a computer cluster with many compute cores, you can use
Message Passing Interface to greatly speed up the modeling process::

    MPI: True
    Nthreads: 1

.. caution:: Nthreads must be equal to 1 if using MPI!

If you want to compare the model results with an image obtained at another
wavelength (e.g., an *HST* image)::

    OpticalImage: /Users/rbussman/Papers/Bussmann_2014a/fitscutouts/XMM03_F110W.fits
    OpticalTag: HST F110W

