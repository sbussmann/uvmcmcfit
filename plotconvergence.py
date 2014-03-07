"""

2013 October 5

Shane Bussmann

Purpose: Plot convergence of lnprob

"""

from astropy.table import Table
import matplotlib.pyplot as plt
from pylab import savefig
import os
import rdconfig
import glob
import numpy

keyname = 'lnprob'

posteriorloc = 'posteriorpdf.dat'
print "Reading burnin results from " + posteriorloc
pdf = Table.read(posteriorloc, format='ascii')
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
