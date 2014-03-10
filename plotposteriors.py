"""
 Plot posterior distributions of model parameters
 Shane Bussmann
 2012 April 5

"""

from __future__ import print_function

#def plot(objname, redo, deltachi2):

import numpy
from astropy.io.misc import hdf5
#import aplpy
import matplotlib.pyplot as mpl
from pylab import savefig
#from matplotlib import rcParams.update({'font.size': 8})
from matplotlib import rc
import modifypdf


rc('font',**{'family':'sans-serif','sans-serif':['Arial Narrow'],'size':'6'})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['New Century Schoolbook']})
#rc('text', usetex=True)

nticks = 5

#deltachi2 = 100
posteriorloc = 'posteriorpdf.hdf5'

# read posterior PDF

print("Reading output from emcee")

fitresults = hdf5.read_table_hdf5(posteriorloc)
fitresults = fitresults[-5000:]
print("prior to pruning <Ln Prob>: {:f}".format(fitresults['lnprob'].mean()))

# identify the good fits
fitresultsgood = modifypdf.prune(fitresults)

# determine dimensions of PDF plots
nparams = len(fitresultsgood[0])
ncol = 4
nrow = nparams / ncol + 1
j = 1

fig = mpl.figure(figsize=(8.0, 1.0 * nrow))

# set up the plotting window
mpl.subplots_adjust(left=0.08, bottom=0.1, right=0.95, top=0.95,
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
            ax = mpl.subplot(nrow, ncol, j)
            j += 1
            mpl.hist(frg, nbins, edgecolor='blue')
            mpl.ylabel('N')
            mpl.xlabel(pname)
            start, end = ax.get_xlim()
            stepsize = (end - start) / nticks
            ax.xaxis.set_ticks(numpy.arange(start, end + 0.99*stepsize, stepsize))

# plot the histogram of chi^2 values
#ax = mpl.subplot(nrow, ncol, j)
#mpl.hist(fitresultsgood['chi2'].max() - fitresultsgood['chi2'], bins=100)
savefig('histposteriors.pdf')
