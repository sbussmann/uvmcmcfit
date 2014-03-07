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

from astropy.io import fits
import numpy
#import pyximport
#pyximport.install(setup_args={"include_dirs":numpy.get_include()})
import sample_vis
import uvutil
#import time


def replace(sbmodelloc, visdataloc):

    # read in the surface brightness map of the model
    modelimage = fits.getdata(sbmodelloc)
    modelheader = fits.getheader(sbmodelloc)

    # read in the uvfits data
    visfile = fits.open(visdataloc)
     
    # load the uv data, including the phase center of the data
    uu, vv = uvutil.uvload(visfile)
     
    # load the uv data, including the phase center of the data
    pcd = uvutil.pcdload(visfile)

    # get the number of visfile[0].data, spws, frequencies, polarizations
    if uu.ndim == 4:
        nvis = uu[:, 0, 0, 0].size
        nspw = uu[0, :, 0, 0].size
        nfreq = uu[0, 0, :, 0].size
        npol = uu[0, 0, 0, :].size
    if uu.ndim == 3:
        nvis = uu[:, 0, 0].size
        nspw = 0
        nfreq = uu[0, :, 0].size
        npol = uu[0, 0, :].size

    # sample the uv visfile[0].data using the model
    uu = uu.flatten()
    vv = vv.flatten()
    model_complex = sample_vis.uvmodel(modelimage, modelheader, \
            uu, vv, pcd)

    if nspw > 0:
        # deflatten the real and imaginary components
        real = numpy.real(model_complex).reshape(nvis, nspw, nfreq, npol)
        imag = numpy.imag(model_complex).reshape(nvis, nspw, nfreq, npol)
        # replace data visfile[0].data with model visfile[0].data
        visfile[0].data['DATA'][:, 0, 0, :, :, :, 0] = real
        visfile[0].data['DATA'][:, 0, 0, :, :, :, 1] = imag
    else:
        # deflatten the real and imaginary components
        real = numpy.real(model_complex).reshape(nvis, nfreq, npol)
        imag = numpy.imag(model_complex).reshape(nvis, nfreq, npol)
        # replace data visfile[0].data with model visfile[0].data
        visfile[0].data['DATA'][:, 0, 0, :, :, 0] = real
        visfile[0].data['DATA'][:, 0, 0, :, :, 1] = imag

    print "Exiting replace"

    #visfile[0].data = vismodel

    return visfile

def subtract(sbmodelloc, visdataloc):

    # read in the surface brightness map of the model
    modelimage = fits.getdata(sbmodelloc)
    modelheader = fits.getheader(sbmodelloc)

    # read in the uvfits data
    visfile = fits.open(visdataloc)
     
    # load the uv data
    uu, vv = uvutil.uvload(visfile)
     
    # load the uv data
    pcd = uvutil.pcdload(visfile)

    # get the number of visfile[0].data, spws, frequencies, polarizations
    if uu.ndim == 4:
        nvis = uu[:, 0, 0, 0].size
        nspw = uu[0, :, 0, 0].size
        nfreq = uu[0, 0, :, 0].size
        npol = uu[0, 0, 0, :].size
    if uu.ndim == 3:
        nvis = uu[:, 0, 0].size
        nspw = 0
        nfreq = uu[0, :, 0].size
        npol = uu[0, 0, :].size

    # sample the uv visfile[0].data using the model
    uu = uu.flatten()
    vv = vv.flatten()
    model_complex = sample_vis.uvmodel(modelimage, modelheader, \
            uu, vv, pcd)

    if nspw > 0:
        # deflatten the real and imaginary components
        real = numpy.real(model_complex).reshape(nvis, nspw, nfreq, npol)
        imag = numpy.imag(model_complex).reshape(nvis, nspw, nfreq, npol)
        # replace data visfile[0].data with model visfile[0].data
        visfile[0].data['DATA'][:, 0, 0, :, :, :, 0] -= real
        visfile[0].data['DATA'][:, 0, 0, :, :, :, 1] -= imag
    else:
        # deflatten the real and imaginary components
        real = numpy.real(model_complex).reshape(nvis, nfreq, npol)
        imag = numpy.imag(model_complex).reshape(nvis, nfreq, npol)
        # replace data visfile[0].data with model visfile[0].data
        visfile[0].data['DATA'][:, 0, 0, :, :, 0] -= real
        visfile[0].data['DATA'][:, 0, 0, :, :, 1] -= imag

    #visfile[0].data = vismodel
    print "Exiting subtract"

    return visfile

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
