Testing Convergence
*******************

The first step in determining whether you have an acceptable model is to plot
the convergence profile of the goodness of fit estimate (by default, the
maximum likelihood estimate).  You can do this using
:func:`visualize.convergence` ::

    import visualize
    visualize.convergence()

This routine produces *convergence.png*, a plot of the difference between
max(lnprob) and lnprob as a function of iteration.  This turns out to be a
convenient way to visualize the convergence, but it does mean that the best-fit
model is shown at the bottom of the plot.

For CDFS_M0, I needed about 450,000 iterations to reach convergence.  The worst
model fit has max(lnprob) - lnprob ~ 2e4.  Both values are fairly typical for a
complex system like CDFS_M0.

.. image::
    media/CDFS_M0_convergence.png
