"""

Created by Shane Bussmann

A number of visualization tools are included here to aid the user in evaluating
the:

    convergence of lnprob
    the posterior PDFs
    the evolution of the PDFs for each parameter of the model
    the covariance matrix for the posterior PDFs
    the best-fit model
    a number of randomly drawn realizations from the posterior PDFs

"""

from __future__ import print_function

import os
from astropy.io import fits
import visualutil
#import sys
#cwd = os.getcwd()
#sys.path.append(cwd)
#import config
import yaml


configloc = 'config.yaml'
configfile = open(configloc)
config = yaml.load(configfile)

def convergence(bestfitloc='posteriorpdf.fits'):

    """

    Plot the convergence profile.  I.e., Max(lnprob) - lnprob as a function of
    iteration number.

    """

    import numpy
    import matplotlib.pyplot as plt
    from pylab import savefig


    print("Reading burnin results from {0:s}".format(bestfitloc))
    pdf = fits.getdata(bestfitloc)
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

def posteriorPDF(bestfitloc='posteriorpdf.fits'):

    """

    Plot the posterior PDF of each parameter of the model.

    """

    # read posterior PDF
    print("Reading output from emcee")
    fitresults = fits.getdata(bestfitloc)
    tag = 'posterior'
    visualutil.plotPDF(fitresults, tag, Ngood=5000, axes='auto')

def evolvePDF(bestfitloc='posteriorpdf.fits', stepsize=50000):

    """

    Plot the evolution of the PDF of each parameter of the model.

    """

    import setuputil


    # Get upper and lower limits on the parameters to set the plot limits
    paramData = setuputil.loadParams(config)
    p_u = paramData['p_u']
    p_l = paramData['p_l']
    limits = [p_l, p_u]

    # read posterior PDF
    fitresults = fits.getdata(bestfitloc)
    nresults = len(fitresults)
    print("Output from emcee has = " + str(nresults) + " iterations.")
    start = 0
    for iresult in range(0, nresults, stepsize):

        strstep = str(stepsize)
        nchar = len(str(nresults))
        striresult = str(iresult).zfill(nchar)
        tag = 'evolution' + strstep + '.' + striresult + '.'
        trimresults = fitresults[start:start + stepsize]
        start += stepsize
        visualutil.plotPDF(trimresults, tag, limits=limits, Ngood=1000, 
                axes='initial')

def covariance(bestfitloc='posteriorpdf.fits'):

    """

    Plot the covariance matrix for the parameters of the model.

    """

    import matplotlib.pyplot as plt
    import numpy
    from pylab import savefig
    import modifypdf
    from astropy.table import Table
    from matplotlib import rc


    # plotting parameters
    rc('font',**{'family':'sans-serif', 'sans-serif':['Arial Narrow'], 
        'size':'6'})

    posteriorpdf = Table.read(bestfitloc)
    posteriorpdf = posteriorpdf[-5000:]

    # remove columns where the values are not changing
    posteriorpdfclean = modifypdf.cleanColumns(posteriorpdf)

    posteriorpdfgood = modifypdf.prune(posteriorpdfclean)

    headers = posteriorpdf.colnames
    ncol = len(headers)
    k = 0
    xsize = ncol * 2
    ysize = ncol * 1.5
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
    savefig('covariance.pdf')
    plt.clf()        

def printFitParam(fitresult, fitKeys, mag=False):
    """ Print parameters for this model                                                         
    mag: bool                                                                                   
        if True, print magnification factors as well                                            
    """
    if mag is False:
        fitresult = fitresult[:-4]
        fitKeys = fitKeys[:-4]

    print("Found the following parameters for this fit:")
    for k, v in zip(fitKeys, fitresult):
        print("%s : %.4f" %(k,v))

def bestFit(bestfitloc='posteriorpdf.fits', showOptical=False, cleanup=True,
        interactive=True):

    """ 

    Read posterior PDF and identify best-fit parameters.  Plot the best-fit
    model and compare to the data.  Also plot the residuals obtained after
    subtracting the best-fit model from the data and compare to the data.
    Optionally plot the best available optical image and compare to the data.
    
    """


    # read the posterior PDFs
    print("Found posterior PDF file: {:s}".format(bestfitloc))
    fitresults = fits.getdata(bestfitloc)

    from astropy.table import Table
    fitKeys = Table.read(bestfitloc).keys()

    # identify best-fit model
    minchi2 = fitresults['lnprob'].max()
    index = fitresults['lnprob'] == minchi2
    bestfit = fitresults[index][0]
    tag = 'bestfit'

    printFitParam(bestfit, fitKeys)
    visualutil.plotFit(config, bestfit, tag=tag, cleanup=cleanup,
            showOptical=showOptical, interactive=interactive)

def goodFits(bestfitloc='posteriorpdf.fits', Nfits=12, Ngood=5000,
        cleanup=True, interactive=True, showOptical=False):

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
    fitresults = fits.getdata(bestfitloc)
    fitresults = fitresults[-Ngood:]
    fitresults = modifypdf.prune(fitresults)
    
    # get keys
    from astropy.table import Table
    fitKeys = Table.read(bestfitloc).keys()

    # select the random realizations model
    Nunprune = len(fitresults)
    realids = numpy.floor(numpy.random.uniform(0, Nunprune, Nfits))

    for ifit in range(Nfits):
        realid = numpy.int(realids[ifit])
        fitresult = fitresults[realid]
        tag = 'goodfit' + str(realid).zfill(4)
        printFitParam(fitresult, fitKeys)
        visualutil.plotFit(config, fitresult, tag=tag, showOptical=showOptical,
                cleanup=cleanup, interactive=interactive)
