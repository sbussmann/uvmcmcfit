"""
2014 February 12

Shane Bussmann

A collection of utilities for modifying posterior probability density
functions.

"""

from astropy.table import Table
import numpy


def trim(oldpdfloc, newpdfloc, niters=5000):

    # get the last niters iterations
    PDFdata = Table.read(oldpdfloc, format='ascii')
    PDFdata = PDFdata[-niters:]

    # write the trimmed list
    PDFdata.write(newpdfloc, format='ascii')

def prune(PDFdata, scaler=5.0):

    # get the last niters iterations
    #PDFdata = Table.read(oldpdfloc, format='ascii')

    print 'prior to pruning: ', PDFdata['lnprob'].mean()
    #import pdb; pdb.set_trace()

    # identify the good fits
    lnprob = PDFdata['lnprob']
    #PDFdata['lnprob'] = lnprob
    minlnprob = lnprob.max()
    dlnprob = numpy.abs(lnprob - minlnprob)
    medlnprob = numpy.median(dlnprob)
    avglnprob = numpy.mean(dlnprob)
    skewlnprob = numpy.abs(avglnprob - medlnprob)
    rmslnprob = numpy.std(dlnprob)
    inliers = (dlnprob < scaler*rmslnprob)
    PDFdata = PDFdata[inliers]
    medlnprob_previous = 0.
    while skewlnprob > 0.1*medlnprob:
        lnprob = PDFdata['lnprob']
        minlnprob = lnprob.max()
        dlnprob = numpy.abs(lnprob - minlnprob)
        rmslnprob = numpy.std(dlnprob)
        inliers = (dlnprob < scaler*rmslnprob)
        PDFdatatmp = PDFdata[inliers]
        if len(PDFdatatmp) == len(PDFdata):
            inliers = (dlnprob < scaler/2.*rmslnprob)
        PDFdata = PDFdata[inliers]
        lnprob = PDFdata['lnprob']
        dlnprob = numpy.abs(lnprob - minlnprob)
        medlnprob = numpy.median(dlnprob)
        avglnprob = numpy.mean(dlnprob)
        skewlnprob = numpy.abs(avglnprob - medlnprob)
        print medlnprob, avglnprob, skewlnprob
        if medlnprob == medlnprob_previous:
            scaler /= 1.5
        medlnprob_previous = medlnprob

    #PDFdatagood = PDFdata.copy()
    # identify the good fits
    goodfits = (PDFdata['lnprob'] <= minlnprob)
    PDFdata = PDFdata[goodfits]

    # return the trimmed list
    #PDFdata.write(newpdfloc, format='ascii')
    return PDFdata
