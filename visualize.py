"""
Created by Shane Bussmann

Purpose: plot 2 things:
1. cleaned map of best-fit model vs. SMA data
2. cleaned map of residual visibilities (model - data)

"""

from __future__ import print_function

import os
from astropy.io.misc import hdf5
import visualutil
import sys
cwd = os.getcwd()
sys.path.append(cwd)
import config


def convergence(bestfitloc='posteriorpdf.hdf5'):

    """

    Plot the convergence profile.  I.e., Max(lnprob) - lnprob as a function of
    iteration number.

    """

    import numpy
    import matplotlib.pyplot as plt
    from pylab import savefig


    print("Reading burnin results from {0:s}".format(bestfitloc))
    pdf = hdf5.read_table_hdf5(bestfitloc)
    keyname = 'lnprob'
    lnprob = pdf[keyname]

    lnprob = numpy.array(lnprob)
    lnprob = lnprob.max() - lnprob
    lnprob = numpy.abs(lnprob)

    plt.clf()
    plt.plot(lnprob, ',', alpha=0.5)
    plt.xlabel('iteration')
    plt.ylabel('max(lnprob) - lnprob')
    tmpcwd = os.getcwd()
    startindx = tmpcwd.find('ModelFits') + 10
    endindx = tmpcwd.find('uvfit') + 7
    objname = tmpcwd[startindx:endindx]
    plt.title(objname)
    plt.semilogy()

    outfile = 'convergence'
    savefig(outfile)

def posteriorPDF(bestfitloc='posteriorpdf.hdf5'):

    """

    Plot the posterior PDF of each parameter of the model.

    """

    import numpy
    import matplotlib.pyplot as plt
    from pylab import savefig
    from matplotlib import rc
    import modifypdf

    rc('font',**{'family':'sans-serif', 'sans-serif':['Arial Narrow'], 
        'size':'6'})

    nticks = 5

    # read posterior PDF

    print("Reading output from emcee")

    fitresults = hdf5.read_table_hdf5(bestfitloc)
    fitresults = fitresults[-5000:]
    print("prior to pruning <Ln Prob>: {:f}".format(fitresults['lnprob'].mean()))

    # identify the good fits
    fitresultsgood = modifypdf.prune(fitresults)

    # determine dimensions of PDF plots
    nparams = len(fitresultsgood[0])
    ncol = 4
    nrow = nparams / ncol + 1
    j = 1

    fig = plt.figure(figsize=(8.0, 1.0 * nrow))

    # set up the plotting window
    plt.subplots_adjust(left=0.08, bottom=0.1, right=0.95, top=0.95,
        wspace=0.4, hspace=0.65)

    pnames = fitresultsgood.keys()

    for pname in pnames:

        # Position of primary lens
        frg = fitresultsgood[pname]
        rmsval = numpy.std(frg)
        if rmsval > 1e-6:
            avgval = numpy.mean(frg)
            print(pname + ' = ', avgval, ' +/- ', rmsval)
            totalwidth = frg.max() - frg.min()
            nbins = totalwidth / rmsval * 5
            ax = plt.subplot(nrow, ncol, j)
            j += 1
            plt.hist(frg, nbins, edgecolor='blue')
            plt.ylabel('N')
            plt.xlabel(pname)
            start, end = ax.get_xlim()
            stepsize = (end - start) / nticks
            ax.xaxis.set_ticks(numpy.arange(start, end + 0.99*stepsize, 
                stepsize))

    savefig('posteriorsPDFs.pdf')

def convariance(bestfitloc='posteriorpdf.hdf5'):

    """

    Plot the covariance matrix for the parameters of the model.

    """

    import matplotlib.pyplot as plt
    import numpy
    from pylab import savefig
    import modifypdf


    posteriorpdf = hdf5.read_table_hdf5(bestfitloc)
    posteriorpdf = posteriorpdf[-5000:]

    posteriorpdfgood = modifypdf.prune(posteriorpdf)

    headers = posteriorpdf.colnames()
    ncol = len(headers)
    k = 0
    xsize = ncol * 3
    ysize = ncol * 2.5
    fig = plt.figure(figsize=(xsize, ysize))
    plt.subplots_adjust(left=0.020, bottom=0.02, right=0.99, top=0.97,
        wspace=0.5, hspace=0.5)

    #for i in numpy.arange(ncol):
    #    ax = plt.subplot(ncol, ncol, i + 1)
    #    namex = 'mu_aper'
    #    namey = headers[i]
    #    datax = mupdfgood[namex]
    #    datay = posteriorpdfgood[namey]
    #    if namex == 'lnprob':
    #        datax = datax.max() - datax
    #    if namey == 'lnprob':
    #        datay = datay.max() - datay
    #    lnprob = posteriorpdfgood['lnprob'].max() - posteriorpdfgood['lnprob']
    #    plt.hexbin(datax, datay, C = lnprob)
    #    plt.xlabel(namex)
    #    plt.ylabel(namey)


    for i in numpy.arange(ncol):

        for j in numpy.arange(ncol - i - 1) + i + 1:

            plotspot = ncol * i + j
            ax = plt.subplot(ncol, ncol, plotspot)

            namex = headers[i]
            namey = headers[j]
            #datax = posteriorpdforig[namex]
            #datay = posteriorpdforig[namey]
            #lnprob = posteriorpdforig['lnprob']

            #plt.hexbin(datax, datay, C = lnprob, color='black')

            datax = posteriorpdfgood[namex]
            datay = posteriorpdfgood[namey]
            if namex == 'lnprob':
                datax = datax.max() - datax
            if namey == 'lnprob':
                datay = datay.max() - datay
            lnprob = posteriorpdfgood['lnprob'].max() - posteriorpdfgood['lnprob']

            plt.hexbin(datax, datay, C = lnprob)
            plt.xlabel(namex)
            plt.ylabel(namey)
            print(i, j, plotspot, namex, namey)
            k += 1

    #plt.suptitle(iau_address, x=0.5, y=0.987, fontsize='xx-large')
    savefig('covariance.png')
    plt.clf()        

def bestFit(bestfitloc='posteriorpdf.hdf5'):

    """ 

    Read posterior PDF and identify best-fit parameters.  Plot the best-fit
    model and compare to the data.  Also plot the residuals obtained after
    subtracting the best-fit model from the data and compare to the data.  
    
    """

    # read the posterior PDFs
    print("Found posterior PDF file: {:s}".format(bestfitloc))
    fitresults = hdf5.read_table_hdf5(bestfitloc)

    # identify best-fit model
    minchi2 = fitresults['lnprob'].max()
    index = fitresults['lnprob'] == minchi2
    bestfit = fitresults[index][0]
    tag = 'bestfit'
    visualutil.preProcess(config, bestfit, tag=tag)

def goodFits(bestfitloc='posteriorpdf.hdf5', Nfits=12, Ngood=5000):

    """ 

    Read posterior PDF and draw Nfits realizations from the final Ngood models
    at random.  Plot the model from each realization and compare to the data.
    Also plot the residuals obtained after subtracting the model from the data
    and compare to the data.  By default: Nfits = 12, Ngood=5000.
    
    """

    import modifypdf
    import numpy


    # read the posterior PDFs
    print("Found posterior PDF file: {:s}".format(bestfitloc))
    fitresults = hdf5.read_table_hdf5(bestfitloc)
    fitresults = fitresults[-Ngood:]
    fitresults = modifypdf.prune(fitresults)

    # select the random realizations model
    Nunprune = len(fitresults)
    realids = numpy.floor(numpy.random.uniform(0, Nunprune, Nfits))

    for ifit in range(Nfits):
        realid = numpy.int(realids[ifit])
        fitresult = fitresults[realid]
        tag = 'goodfit' + str(realid).zfill(4)
        visualutil.preProcess(config, fitresult, tag=tag)
