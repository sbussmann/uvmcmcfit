"""

Inspect a model with specific parameters defined in config_sandbox.py

"""


from __future__ import print_function

import os
import visualutil
import sys
cwd = os.getcwd()
sys.path.append(cwd)
import config_sandbox
import setuputil
import numpy


def plot():

    # read the input parameters
    paramData = setuputil.loadSandboxParams(config_sandbox)
    parameters_0 = numpy.array(paramData['parameters'])
    nregions = paramData['nregions']

    # Loop over each region
    npar_previous = 0
    for regioni in range(nregions):

        # search poff_models for parameters fixed relative to other parameters
        fixindx = setuputil.fixParams(paramData)
        poff = paramData['poff']
        ndim_total = len(poff)
        fixed = (numpy.where(fixindx >= 0))[0]
        nfixed = fixindx[fixed].size
        parameters_offset = numpy.zeros(ndim_total)
        for ifix in range(nfixed):
            ifixed = fixed[ifix]
            subindx = fixindx[ifixed]
            par0 = 0
            if fixindx[subindx] > 0:
                par0 = parameters_0[fixindx[subindx]]
            parameters_offset[ifixed] = parameters_0[subindx] + par0

        parameters_regions = parameters_0 + parameters_offset

        nlens = config_sandbox.Nlens[regioni]
        nsource = config_sandbox.Nsource[regioni]

        nparlens = 5 * nlens
        nparsource = 6 * nsource
        npar = nparlens + nparsource + npar_previous
        parameters = parameters_regions[npar_previous:npar]
        npar_previous = npar

        visualutil.plotFit(config_sandbox, paramData, parameters, regioni, 
                tag='sandbox', cleanup=True)
