"""

2015 March 24

Shane Bussmann

Compare simulated visibilities from my implementation of miriad's uvmodel
routine and to miriad's implementation.  Also plot the difference.  Do this for
the real and imaginary components as well as the weights.  And possibly the u,
v, and w values as well.

"""

import uvutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import savefig
from subprocess import call
import uvmodel


def iterPlot(simvis, nfigures, color, label):
    """
    Plot simulated visibilities, weights, and uvw values vs. entry number.
    """

    for i in range(nfigures):
        plt.figure(i)
        plt.plot(simvis[i, :], ',', color=color, label=label)

def rollUp(visfile):
    u, v, w = uvutil.uvload(visfile)
    data, weight = uvutil.visload(visfile)
    real = np.real(data)
    imag = np.imag(data)
    nvis = u.size
    rolled = np.zeros((6, nvis))
    rolled[0, :] = u.flatten()
    rolled[1, :] = v.flatten()
    rolled[2, :] = w.flatten()
    rolled[3, :] = real.flatten()
    rolled[4, :] = imag.flatten()
    rolled[5, :] = weight.flatten()
    return rolled

def savePlots(names, nfigures, colors):
    for i in range(nfigures):
        name = names[i]
        plt.figure(i)
        plt.ylabel(name)
        plt.xlabel('Observation Number')
        plt.legend()
        leg = plt.legend(loc = 'best')
        for i, text in enumerate(leg.get_texts()):
            color = colors[i]
            plt.setp(text, color = color)
        plt.tight_layout()
        savefig('viscompare_' + name + '.pdf')
        print("Finished saving " + name + " figure!")

def iterFig(uvmcmcfitFile, miriadFile):

    uvmcmcfit = rollUp(uvmcmcfitFile)
    miriad = rollUp(miriadFile)
    #miriad[5, :] = miriad[5, :] * 2
    difference = uvmcmcfit - miriad

    nfigures = uvmcmcfit[:, 0].size

    plt.ioff()
    iterPlot(uvmcmcfit, nfigures, 'red', 'uvmcmcfit')
    print("Finished plotting uvmcmcfit data!")
    iterPlot(miriad, nfigures, 'blue', 'miriad')
    print("Finished plotting miriad data!")
    iterPlot(difference, nfigures, 'black', 'difference')
    print("Finished plotting uvmcmcfit - miriad data!")

    names = ['u', 'v', 'w', 'real', 'imaginary', 'weights']
    colors = ['red', 'blue', 'black']
    savePlots(names, nfigures, colors)

def miriadVis(model, data, simfile):
    """
    Use miriad to make a simulated visibility dataset given a model and data.
    """

    # first turn the model image into miriad format
    try:
        index = model.index('fits')
        modelmiriad = model[0:index] + 'miriad'
        cmd = 'rm -rf ' + modelmiriad
        call(cmd, shell=True)
    except:
        print("Uh oh, couldn't remove the existing file.")
    try:
        cmd = 'fits op=xyin in=' + model + ' out=' + modelmiriad
        call(cmd, shell=True)
    except:
        print("Uh oh, couldn't make the miriad model image.")

    # next, turn the observed visibilities into miriad format
    index = data.index('uvfits')
    datamiriad = data[0:index] + 'miriad'
    cmd = 'rm -rf ' + datamiriad
    call(cmd, shell=True)
    cmd = 'fits op=uvin in=' + data + ' out=' + datamiriad
    call(cmd, shell=True)

    # next, make the simulated visibilities
    cmd = 'rm -rf ' + simfile
    call(cmd, shell=True)
    cmd = 'uvmodel options=replace model=' + modelmiriad + ' vis=' + \
            datamiriad + ' out=' + simfile
    call(cmd, shell=True)

    # and convert the simulated visibilities to uvfits format
    index = simfile.index('miriad')
    simfilefits = simfile[0:index] + 'uvfits'
    cmd = 'rm -rf ' + simfilefits
    call(cmd, shell=True)
    cmd = 'fits op=uvout in=' + simfile + ' out=' + simfilefits
    call(cmd, shell=True)

def uvmcmcfitVis(model, data, simfile):
    #cmd = 'rm -rf ' + simfile
    #call(cmd, shell=True)
    uvmodel.replace(model, data, simfile, miriad=True)
