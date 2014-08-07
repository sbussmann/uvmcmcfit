Accessing the Data Directly
***************************

The model fit results are stored in a fits table called *posteriorpdf.fits*.
You can inspect the results directly using the :mod:`astropy.table` module::

    # import astropy's table module
    import astropy.table as Table

    # read the fit results
    fitresults = Table.read('posteriorpdf.fits')

    # get and print the column header names
    keys = fitresults.keys()
    print(keys)

    # you can plot whatever aspect of fitresults you want from here...
