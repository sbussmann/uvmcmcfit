"""

2013 August 8
Shane Bussmann

Purpose: Look for degeneracies between parameters of lens model

"""

import matplotlib.pyplot as plt
from astropy.table import Table
import numpy
from pylab import savefig
import modifypdf


# examine each inputs_lens.dat file to determine total number of lenses

posteriorpdfloc = 'posteriorpdf.dat'
posteriorpdf = Table.read(posteriorpdfloc, format='ascii', data_start=-5000)

posteriorpdfgood = modifypdf.prune(posteriorpdf)

headers = posteriorpdf.keys()
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
        print i, j, plotspot, namex, namey
        k += 1

#plt.suptitle(iau_address, x=0.5, y=0.987, fontsize='xx-large')
savefig('degeneracies.png')
plt.clf()

