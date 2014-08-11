Plot the Covariance Matrix
**************************

You might wish to examine how the various parameters are correlated with each
other.  You can do this by plotting the covariance matrix using
:func:`visualize.covariance`::

    import visualize
    visualize.covariance()

This will produce an enormous plot showing correlations between all of the
model parameters.  The plot is too big to reproduce here, but you can use
Preview to zoom in on any panel of interest.  There is evidence that
``DeltaRA`` for ``Lens0`` is highly anti-correlated with ``DeltaRA`` for
``Source0``.  This is expected since the source position has been fixed
relative to the lens position in this model fit.  The same holds true for
``DeltaDec``.
