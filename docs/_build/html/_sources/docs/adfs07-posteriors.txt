Plot the Posterior PDFs
***********************

You can examine the posterior probability density functions to get a
quantitative look at the results of the model fitting process using
:func:`visualize.posteriorPDF`::

    import visualize
    visualize.posteriorPDF()

This will produce a series of histograms showing the posterior probability
distribution functions for every parameter in the model.

.. image::
    ADFS07_posteriorPDFs.png

This routine also prints the average and 1-sigma rms uncertainty on each
parameter of the model.  

You can also see how the posterior PDF of every parameter in the
model changes as a function of iteration using :func:`visualize.evolvePDF`::

    visualize.evolvePDF()

This function essentially produces a posteriorPDF every ``stepsize``
iterations.  The default is ``stepsize = 50000``.  You can then view the
evolution in the PDF using a viewer application like Preview.
