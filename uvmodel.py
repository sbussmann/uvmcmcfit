"""
2013 May 13
Shane Bussmann

Purpose: Take an input uv fits file and replace or subtract the visibilities
with a model for the surface brightness map of the target.

Inputs: 
    sbmodelloc: location of fits image of model for surface brightness map
    visdataloc: location of uv fits data for target

Outputs:
    vismodel: uvfits data for model
    visheader: uvfits header for model
"""

from __future__ import print_function

from astropy.io import fits
import numpy
#import pyximport
#pyximport.install(setup_args={"include_dirs":numpy.get_include()})
import sample_vis
import uvutil
import os
#import time


def writeVis(vis_complex, visdataloc, modelvisloc, miriad=False):

    if miriad:
        os.system('rm -rf ' + modelvisloc)
        cmd = 'cp ' + visdataloc + ' ' + modelvisloc
        os.system(cmd)
        # get the real and imaginary components
        real = numpy.real(vis_complex)
        imag = numpy.imag(vis_complex)

        # replace data visibilities with model visibilities
        visfile = fits.open(modelvisloc, mode='update')
        visibilities = visfile[0].data

        visheader = visfile[0].header
        if visheader['NAXIS'] == 7:
            visibilities['DATA'][:, 0, 0, :, :, :, 0] = real
            visibilities['DATA'][:, 0, 0, :, :, :, 1] = imag
        elif visheader['NAXIS'] == 6:
            visibilities['DATA'][:, 0, 0, :, :, 0] = real
            visibilities['DATA'][:, 0, 0, :, :, 1] = imag
        else:
            print("Visibility dataset has >7 or <6 axes.  I can't read this.")
        #visibilities['DATA'][:, 0, 0, :, :, 2] = wgt

        # replace the data visibilities with the model visibilities
        visfile[0].data = visibilities
        #visfile.flush()
        
    else:
        from taskinit import tb
        print(modelvisloc)
        tb.open(visdataloc)
        os.system('rm -rf ' + modelvisloc)
        tb.copy(modelvisloc)
        tb.close()
        tb.open(modelvisloc, nomodify=False)
        datashape = tb.getcol('DATA').shape
        vis_complex = vis_complex.reshape(datashape)
        tb.putcol('DATA', vis_complex)
        tb.close()

def getVis(sbmodelloc, visdataloc):

    # read in the surface brightness map of the model
    modelimage = fits.getdata(sbmodelloc)
    modelheader = fits.getheader(sbmodelloc)
     
    # load the uv data, including the phase center of the data
    uu, vv = uvutil.uvload(visdataloc)
     
    # load the uv data, including the phase center of the data
    pcd = uvutil.pcdload(visdataloc)

    # sample the uv visfile[0].data using the model
    uushape = uu.shape
    if len(uushape) == 2:
        npol = uushape[0]
        nrow = uushape[1]
        uushape = (npol, 1, nrow)
    uu = uu.flatten()
    vv = vv.flatten()
    model_complex = sample_vis.uvmodel(modelimage, modelheader, \
            uu, vv, pcd)

    vis_complex = model_complex.reshape(uushape)
    return vis_complex

def replace(sbmodelloc, visdataloc, modelvisloc, miriad=False):

    vis_model = getVis(sbmodelloc, visdataloc)
    #sub_complex, vis_weight = uvutil.visload(visdataloc)

    writeVis(vis_model, visdataloc, modelvisloc, miriad=miriad)

    return

def subtract(sbmodelloc, visdataloc, modelvisloc, miriad=False):

    vis_model = getVis(sbmodelloc, visdataloc)
     
    # load the visibilities
    vis_data, vis_weight = uvutil.visload(visdataloc)

    vis_data -= vis_model

    #print(visdataloc, modelvisloc, miriad)
    writeVis(vis_data, visdataloc, modelvisloc, miriad=miriad)

    return

def add(sbmodelloc, visdataloc, WeightByRMS=True, ExcludeChannels='none'):

    # read in the surface brightness map of the model
    modelimage = fits.getdata(sbmodelloc)
    modelheader = fits.getheader(sbmodelloc)

    # read in the uvfits data
    visfile = fits.open(visdataloc)
    visibilities = visfile[0].data
     
    # load the uv data, including the phase center of the data
    uu, vv, pcd = uvutil.uvload(visfile)

    # get the number of visibilities, spws, frequencies, polarizations
    if uu.ndim == 4:
        nvis = uu[:, 0, 0, 0]
        nspw = uu[0, :, 0, 0]
        nfreq = uu[0, 0, :, 0]
        npol = uu[0, 0, 0, :]
    if uu.ndim == 3:
        nvis = uu[:, 0, 0]
        nspw = 0
        nfreq = uu[0, :, 0]
        npol = uu[0, 0, :]

    # sample the uv visibilities using the model
    uu = uu.flatten()
    vv = vv.flatten()
    model_complex = sample_vis.uvmodel(modelimage, modelheader, \
            uu, vv, pcd)

    if nspw > 0:
        # get the real and imaginary components
        real = numpy.real(model_complex).reshape(nvis, nspw, nfreq, npol)
        imag = numpy.imag(model_complex).reshape(nvis, nspw, nfreq, npol)
        uu = uu.reshape(nvis, nspw, nfreq, npol)
        vv = vv.reshape(nvis, nspw, nfreq, npol)

        # replace data visibilities with model visibilities
        visibilities['DATA'][:, 0, 0, :, :, :, 0] = real_raw + real
        visibilities['DATA'][:, 0, 0, :, :, :, 1] = imag_raw + imag
        visibilities['DATA'][:, 0, 0, :, :, :, 2] = wgt
    else:
        # get the real and imaginary components
        real = numpy.real(model_complex).reshape(nvis, nfreq, npol)
        imag = numpy.imag(model_complex).reshape(nvis, nfreq, npol)
        uu = uu.reshape(nvis, nfreq, npol)
        vv = vv.reshape(nvis, nfreq, npol)

        # replace data visibilities with model visibilities
        visibilities['DATA'][:, 0, 0, :, :, 0] = real_raw + real
        visibilities['DATA'][:, 0, 0, :, :, 1] = imag_raw + imag
        visibilities['DATA'][:, 0, 0, :, :, 2] = wgt

    # replace the data visibilities with the model visibilities
    visfile[0].data = visibilities

    #visfile[0].data = vismodel

    return visfile
