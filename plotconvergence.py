"""

2013 October 5

Shane Bussmann

Purpose: Plot convergence of lnprob

"""

from astropy.io.misc import hdf5
import matplotlib.pyplot as plt
from pylab import savefig
import os
import numpy

keyname = 'lnprob'

posteriorloc = 'posteriorpdf.hdf5'
print "Reading burnin results from " + posteriorloc
pdf = hdf5.read_table_hdf5(posteriorloc)
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
