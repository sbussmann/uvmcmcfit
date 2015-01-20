"""

Inspect a model with specific parameters defined in config_sandbox.py

"""


from __future__ import print_function

import visualutil
import setuputil
import yaml
import numpy


def plot(cleanup=True, configloc='sandbox.yaml', interactive=True):

    # read the input parameters
    configfile = open(configloc, 'r')
    config = yaml.load(configfile)
    
    paramData = setuputil.loadParams(config)
    testfit = paramData['p_l']
    # prepend 1 dummy value to represent lnprob
    testfit = numpy.append(0, testfit)
    # append nmu dummy values that represent the magnification factors
    nlensedsource = paramData['nlensedsource']
    nlensedregions = paramData['nlensedregions']
    nmu = 2 * (numpy.array(nlensedsource).sum() + nlensedregions)
    
    for i in range(nmu):
        testfit = numpy.append(testfit, 0)
    tag = 'sandbox'
    visualutil.plotFit(config, testfit, tag=tag, cleanup=cleanup,
            interactive=interactive)
        
