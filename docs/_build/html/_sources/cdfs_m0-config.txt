Configuring config.yaml
***********************

config.yaml contains the instructions needed by ``uvmcmcfit`` to initiate the
model fitting process.

Required keywords
-----------------

A few house-keeping parameters::

    # Name of the target
    ObjectName: CDFS_M0
 
    # Name of the fits image; the pixel scale in this image sets the pixel
    # scale in the model image
    ImageName: CDFS_M0.Xabf.selfcal.statwt.cont.mfs.fits

    # Name of the uvfits visibility data; the weights should be scaled such
    # that Sum(weights * real) ~ N_vis [see uvutil.statwt()]
    UVData: CDFS_M0.Xabf.selfcal.statwt.cont.uvfits
    
    # Number of walkers
    Nwalkers: 92

.. caution:: 

    The number of walkers used by emcee must be more than double the number of
    parameters).  In this case, there are 23 parameters (5 for the lens, 6 for
    each of three sources), so the minimum number of walkers is 46.  I selected
    92 to be on the safe side.

Now for parameters that describe the geometry of the system.  You must define
at least one region.  The first region should be named ``Region0``, the second
``Region1``, etc.  Pay attention to the indentation; the remaining keywords
must be indented to indicate they are sub-components of ``Region0``.  For each
region, you must define a RA and Dec center, an angular radial extent that
contains the emission which you are attempting to model, and at least one
source.  You have the option to use the ``Oversample`` keyword to specify an
integer factor by which to increase the resolution of the model image relative
to the data image (i.e., relative to the resolution in the image specified by
``ImageName``).  

The first source should be named ``Source0``, the second source should be named
``Source1``, etc.  Sources are elliptical Gaussians.  Each source must have the
following parameters: the total intrinsic flux density (IntrinsicFlux [mJy]),
the effective radius (EffectiveRadius, defined as sqrt(a*b) [arcsec]), the
offset in RA and Dec from RACentroid and DecCentroid (DeltaRA and DeltaDec
[arcsec]), the axial ratio (AxialRatio, defined as b/a), and the position angle
in degrees east of north (PositionAngle [degrees]).

For each source parameter, you must specify the lower and upper limits as well
as how to initialize the walkers for that parameter.  This is done using the
following syntax: ``Limits: [lower limit, lower initialization, upper
initialization, upper limit]``. So, for example, in the code snippet below for
CDFS_M0, ``Source0`` is permitted to have a total intrinsic flux density ranging
from 0.1 to 30 mJy, but is initialized with a uniform probability density
distribution between 1 and 5 mJy.

You may account for the deflection of light rays caused by the presence of a
galaxy or group of galaxies acting as a gravitational lens by specifying one or
more lenses.  The first lens should be named ``Lens0``, the second lens should
be named ``Lens1``, etc.  

Lenses are assumed to be singular isothermal ellipsoids.  They are
parameterized by: the Einstein radius (EinsteinRadius [arcsec]), the offset in
RA and Dec from RACentroid and DecCentroid (DeltaRA and DeltaDec [arcsec]), the
axial ratio (AxialRatio), and the position angle in degrees east of north
(PositionAngle [degrees]).

Lens parameters are specified in the same way as source parameters: ``Limits:
[lower limit, lower initialization, upper initialization, upper limit]``.  

.. Note:: 

    It is sometimes desirable to specify the permissible range on a given
    parameter relative to another parameter of the model.  For example, you
    might wish to force ``Source0`` to be north of ``Lens0``.  You can
    accomplish this by adding a line under the ``Limits`` specification for
    ``DeltaDec`` for ``Source0`` that says ``FixedTo: Region0 Lens0 DeltaDec``.
    This makes the program understand that DecSource0 = DecCentroid +
    DeltaDecLens0 + DeltaDecSource0, rather than simply DecSource0 =
    DecCentroid + DeltaDecSource0.  The example below shows how to fix the RA
    and Dec of ``Source0`` relative to the RA and Dec of ``Lens0``.

.. Note:: 

    In some cases, both the lens itself and the lensed emission is detected by
    the interferometer.  The best way I have found to deal with this situation
    is to create two regions with the same position center and angular extent,
    ``Region0`` and ``Region1``, corresponding to the lensed and lens emission,
    respectively.  These regions will be modeled simultaneously, so that there
    is no need to do a "lens-subtraction" prior to modeling the lensed
    emission.  An example based on CDFS_M0 follows below.

.. Note:: 

    You can fix a parameter in the model to a given value by specifying both
    the lower and upper initialization to have that same value.  An example of
    how to do this is shown for ``DeltaRA`` and ``DeltaDec`` in ``Source0`` of
    ``Region1``.

::

    # First region: contains emission from the lensed galaxies
    Region0:

        # Right Ascension and Declination center of the model image (degrees)::
        RACentroid: 51.966712
        DecCentroid: -29.152889

        # Angular radial extent of the model image (arcsec)
        RadialExtent: 3.0

        # [OPTIONAL]
        # Integer factor by which to increase resolution of model image
        Oversample: 2

        # Source0 -- this source is strongly lensed
        Source0:

            # total intrinsic flux density
            IntrinsicFlux:
                Limits: [0.1, 1.0, 2.0, 30.0]

            # effective radius of elliptical Gaussian [sqrt(a*b)] (arcsec)
            EffectiveRadius:
                Limits: [0.01, 0.1, 0.4, 1.5]

            # Offset in RA and Dec from RALens0 and DecLens0 (arcseconds)
            DeltaRA:
                FixedTo: Region0 Lens0 DeltaRA
                Limits: [-1.7, -0.4, -0.3, 1.7]
            DeltaDec:
                FixedTo: Region0 Lens0 DeltaDec
                Limits: [-0.7, 0.1, 0.2, 0.7]

            # axial ratio = semi-minor axis / semi-major axis
            AxialRatio:
                Limits: [0.3, 0.3, 1.0, 1.0]

            # position angle (degrees east of north)
            PositionAngle:
                Limits: [0.0, 0.0, 180.0, 180.0]

        # Source1 -- this source is weakly lensed
        Source1:

            # total intrinsic flux density
            IntrinsicFlux:
                Limits: [0.1, 6.0, 8.0, 30.0]

            # effective radius of elliptical Gaussian [sqrt(a*b)] (arcsec)
            EffectiveRadius:
                Limits: [0.01, 0.1, 0.4, 1.5]

            # Offset in RA and Dec from RACentroid and DecCentroid (arcseconds)
            DeltaRA:
                Limits: [-1.2, -0.5, 0.0, 0.2]
            DeltaDec:
                Limits: [0.5, 1.2, 1.8, 2.5]

            # axial ratio = semi-minor axis / semi-major axis = b/a
            AxialRatio:
                Limits: [0.3, 0.3, 1.0, 1.0]

            # position angle (degrees east of north)
            PositionAngle:
                Limits: [0.0, 0.0, 180.0, 180.0]

        # Lens0
        Lens0:

            # Einstein radius
            EinsteinRadius:
                Limits: [0.4, 1.0, 1.5, 2.0]

            # Offset in RA and Dec from RACentroid and DecCentroid (arcseconds)
            DeltaRA:
                Limits: [0.1, 0.2, 0.25, 0.3]
            DeltaDec:
                Limits: [-1.9, -1.8, -1.75, -1.7]

            # axial ratio = semi-minor axis / semi-major axis
            AxialRatio:
                Limits: [0.3, 0.7, 0.9, 1.0]

            # position angle (degrees east of north)
            PositionAngle:
                Limits: [0.0, 0.0, 180.0, 180.0]

    # Second region: contains emission from the lens
    Region1:

        # Right Ascension and Declination center of the model image (degrees)::
        RACentroid: 51.966712
        DecCentroid: -29.152889

        # Angular radial extent of the model image (arcsec)
        RadialExtent: 3.0

        # [OPTIONAL]
        # Integer factor by which to increase resolution of model image
        Oversample: 2

        # Source0 -- this is the lens
        Source0:

            # total intrinsic flux density
            IntrinsicFlux:
                Limits: [0.1, 5.0, 6.0, 30.0]

            # effective radius of elliptical Gaussian [sqrt(a*b)] (arcsec)
            EffectiveRadius:
                Limits: [0.01, 0.1, 0.2, 0.5]

            # Offset in RA and Dec from RALens0 and DecLens0 (arcseconds)
            # I assume the center of the gravitational potential is coincident
            # with the emission centroid of the lensing galaxy
            DeltaRA:
                FixedTo: Region0 Lens0 DeltaRA
                Limits: [-1.7, 0.0, 0.0, 1.7]
            DeltaDec:
                FixedTo: Region0 Lens0 DeltaDec
                Limits: [-0.7, 0.0, 0.0, 0.7]

            # axial ratio = semi-minor axis / semi-major axis
            AxialRatio:
                Limits: [0.3, 0.3, 1.0, 1.0]

            # position angle (degrees east of north)
            PositionAngle:
                Limits: [0.0, 0.0, 180.0, 180.0]

      

Optional keywords
-----------------

By default, the maximum likelihood estimate is used to measure the goodness of
fit.  Alternatively, you may use the chi-squared value as the goodness of fit
criterion via::

    # Goodness of fit measurement
    LogLike: chi2

By default, parallel processing is not used.  To use parallel processing on a
single machine, set the Nthreads variable to a number greater than 1.  For
example, ::

    # Number of threads for multi-processing on a single computer
    Nthreads: 2

If you have access to a computer cluster with many compute cores, you can use
Message Passing Interface to greatly speed up the modeling process::

    # Use Message Passing Interface
    MPI: True
    Nthreads: 1

.. caution:: Nthreads must be equal to 1 if using MPI!

If you want to compare the model results with an image obtained at another
wavelength (e.g., an *HST* image), you must specify the location of the
alternative image as well as the telescope and filter used to obtain the
image::

    # Alternative image name (used only for comparing with best-fit model)
    OpticalImage: CDFS_M0_Ks.fits

    # Telescope and filter of alternative image    
    OpticalTag: VIDEO Ks

