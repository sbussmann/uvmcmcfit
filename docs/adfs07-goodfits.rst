Plot Several Acceptable Model Fits
**********************************

It's nice to have visual confirmation that the best-fit model gives an
acceptable fit to the data.  It's even better to see that a random draw from
the posterior probability density function (PDF) also gives an acceptable fit
to the data.  This is easily done using::

    import visualize
    visualize.goodFits()

By default, :func:`visualize.goodFits` produces plots of 12 model fits randomly
drawn from the posterior PDF.  You can adjust this with the ``Nfits`` keyword
argument.

.. note::

    Other than ``Nfits``, :func:`visualize.goodFits` takes the same optional
    keyword arguments as :func:`visualize.bestFit`.

