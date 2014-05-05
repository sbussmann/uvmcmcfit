"""
Created by Shane Bussmann

Purpose: plot 2 things:
1. cleaned map of best-fit model vs. SMA data
2. cleaned map of residual visibilities (model - data)

"""

#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['New Century Schoolbook']})
#rc('text', usetex=True)


import math
import numpy
import os
from astropy import wcs
from astropy.io import fits
from astropy.io.misc import hdf5
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from pylab import savefig
import lensutil
import uvmodel
import uvutil
import sys
cwd = os.getcwd()
sys.path.append(cwd)
import config
#from scipy import ndimage
#import matplotlib.font_manager as fm


# set font properties
font = {'family' : 'Arial Narrow',
        'weight' : 'bold',
        'size'   : 10}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 1.5

# get the current working directory
cwd = os.getcwd()

# Specify the location of the output files
outfits_miriad = 'sbmap_image'
outfits_aper = 'aperture_image.fits'
outfits_unlens = 'sbmap_source'
outfits_unlens_aper = 'aperture_source.fits'
outdef1 = 'outdef1.dat'
outcrit = 'outcrit.dat'
#vis_mir = 'sbmap_vis.fits'

yes_shear = False

#--------------------------------------------------------------------------
# Step 1: read in interferometric image and beam
im = fits.getdata(config.ImageName)
im = im[0,0,:,:].copy()
headim = fits.getheader(config.ImageName)

# get the object name from directory
#startindx = cwd.find('ModelFits')
#endindx = cwd.find('uvfit')
#objectname = cwd[startindx + 10 : endindx - 1]
try:
    objectname = headim['OBJECT']
except KeyError:
    objectname = 'No objname in header'

# Obtain measurements of beamsize and image min/max
bmaj = headim['BMAJ'] * 3600
bmin = headim['BMIN'] * 3600
#bmaj2 = headim2['BMAJ'] * 3600
#bmin2 = headim2['BMIN'] * 3600
bpa = headim['BPA']
cdelt1 = headim['CDELT1'] * 3600
cdelt2 = headim['CDELT2'] * 3600
celldata = math.sqrt( abs(cdelt1) * abs(cdelt2) )
cdelt = celldata / 3600 #* math.pi / 180
strcdelt = str(cdelt)
npix = math.pi * bmaj/2 * bmin/2 / celldata**2 / math.log(2)

datawcs = wcs.WCS(headim, naxis=2)
nx = im[0,:].size
ny = im[:,0].size
xcenter = nx / 2 + 1
ycenter = ny / 2 + 1
#sky = wcs.wcs_pix2sky(xcenter, ycenter, 1)
#ra_center = sky[0]
#dec_center = sky[1]

# compute rms within central 3/4 of the image
mask = im.copy()
mask[:] = 1
yr0 = ny / 4
yr1 = 3 * ny / 4
xr0 = nx / 4
xr1 = 3 * nx / 4
mask[yr0:yr1, xr0:xr1] = 0

# determine the number of regions for which we need surface brightness maps
regionIDs = config.RegionID
nregions = len(regionIDs)

# determine the number of models which are lensed
#contents = rdconfig.LineContents(configloc, 'Lensed')
#arraycontents = numpy.array(contents)
#lensed = (numpy.where(arraycontents == 'True'))[0]
#nlensmodel = lensed.size

nlens_regions = []
nsource_regions = []
poff = []
pname = []
model_types_regions = []
for regioni in range(nregions):
    ri = str(regioni)
    ra_centroid = config.RACentroid[regioni]
    dec_centroid = config.DecCentroid[regioni]
    extent = config.RadialExtent[regioni]
    oversample = config.Oversample[regioni]
    nlens = config.Nlens[regioni]
    nsource = config.Nsource[regioni]

    # Append the number of lenses and sources for this region
    nlens_regions.append(nlens)
    nsource_regions.append(nsource)

    for ilens in range(nlens):

        li = str(ilens)
        lensparams = ['EinsteinRadius', 'DeltaRA', 'DeltaDec', 'AxialRatio', \
                'PositionAngle']
        tag = '_Lens' + li + '_Region' + ri
        for lensparam in lensparams:
            fullparname = 'Constraint_' + lensparam + tag
            values = getattr(config, fullparname)
            poff.append(values.pop()) 
            values = numpy.array(values).astype(float)
            pname.append(lensparam + tag)
            fullparname = 'Init_' + lensparam + tag
            values = getattr(config, fullparname)

    for isource in range(nsource):

        si = str(isource)

        sourceparams = ['IntrinsicFlux', 'Size', 'DeltaRA', 'DeltaDec', \
                'AxialRatio', 'PositionAngle']
        tag = '_Source' + si + '_Region' + ri
        for sourceparam in sourceparams:
            fullparname = 'Constraint_' + sourceparam + tag
            values = getattr(config, fullparname)
            poff.append(values.pop()) 
            values = numpy.array(values).astype(float)
            pname.append(sourceparam + tag)
            fullparname = 'Init_' + sourceparam + tag
            values = getattr(config, fullparname)

        # get the model type
        fullparname = 'ModelMorphology' + tag
        model_types_regions.append(getattr(config, fullparname))

    # mask regions containing significant emission
    skyxy = datawcs.wcs_world2pix(ra_centroid, dec_centroid, 1)
    x_center = skyxy[0]
    y_center = skyxy[1]
    pixextent = extent / celldata
    xm0 = x_center - pixextent / 2
    xm1 = x_center + pixextent / 2
    ym0 = y_center - pixextent / 2
    ym1 = y_center + pixextent / 2
    mask[ym0:ym1, xm0:xm1] = 2

# determine the indices for fixed parameters
ndim_total = len(poff)
fixindx = numpy.zeros(ndim_total) - 1
for ifix in range(ndim_total):
    if pname.count(poff[ifix]) > 0:
        fixindx[ifix] = pname.index(poff[ifix])

# good region is where mask is zero
goodregion = mask == 0

# compute sigma image from cutout of SMA flux image
#Y = histogram( im1[goodregion], bin=0.0006, locations=X )
#Result = GAUSSFIT( X, Y, A )
#rms = A[2]
rms = im[goodregion].std()
#print rms
#npix_sma2 = math.pi * bmaj2/2 * bmin2/2 / celldata**2 / math.log(2)
immin = -rms
immax = im.max()

#------------------------------------------------------------------------------
# Read best-fit results file
bestfitloc = 'posteriorpdf.hdf5'

# read the latest posterior PDFs
print "Found latest posterior PDF file: " + bestfitloc
fitresults = hdf5.read_table_hdf5(bestfitloc)
#fitresults = Table.read(bestfitloc, format='ascii')

# identify best-fit model
minchi2 = fitresults['lnprob'].max()
indx = fitresults['lnprob'] == minchi2#fitresults['lnprob'][1]
bestfit = fitresults[indx][0]

nmu = 2 * (numpy.array(nsource_regions).sum() + nregions)
pzero_regions = list(bestfit.data)[1:-nmu]
#if len(bestfit) > 1:
#bestfit = bestfit[0]

print objectname, bestfit.data
#rint bestfit['shear']

# search poff_models for parameters fixed relative to other parameters
fixed = (numpy.where(fixindx >= 0))[0]
nfixed = fixindx[fixed].size
poff_regions = numpy.zeros(ndim_total)
for ifix in range(nfixed):
    poff_regions[fixed[ifix]] = bestfit[fixindx[fixed[ifix]] + 1]

parameters_regions = pzero_regions + poff_regions

#-----------------------------------------------------------------
# Skip gravlens:
# Create a surface brightness map of lensed emission for the given set
# of foreground lens(es) and background source parameters.
#-----------------------------------------------------------------

# hack to get the right index from the pzero variable
ind = 0
npar_previous = 0
prindx = 0

# loop over each model
for regioni in range(nregions):
    ri = str(regioni)
    sri = str(ri)
    ra_centroid = config.RACentroid[regioni]
    dec_centroid = config.DecCentroid[regioni]
    extent = config.RadialExtent[regioni]
    oversample = config.Oversample[regioni]
    nlens = config.Nlens[regioni]
    nsource = config.Nsource[regioni]
    model_types = model_types_regions[prindx:prindx + nsource]

    nparlens = 5 * nlens
    nparsource = 6 * nsource
    npar = nparlens + nparsource + npar_previous
    parameters = parameters_regions[npar_previous:npar]
    npar_previous += npar

    # define spatial range over which lensed emission is detected
    xlo = -extent
    xhi = extent
    ylo = -extent
    yhi = extent

    # define number of pixels in lensed surface brightness map
    dx = xhi - xlo
    nxmod = oversample * int( round( dx / celldata ) )
    dy = yhi - ylo
    nymod = oversample * int( round( dy / celldata ) )

    # make x and y coordinate images for lens model
    x = dx * numpy.outer(numpy.ones(nymod), numpy.arange(nxmod)) / \
            float(nxmod - 1) + xlo
    y = dy * numpy.outer(numpy.arange(nymod), numpy.ones(nxmod)) / \
            float(nymod - 1) + ylo

    # Provide world-coordinate system transformation data in the header of
    # the lensed surface brightness map
    modelheader = headim.copy()
    crpix1 = nxmod / 2 + 1
    crpix2 = nymod / 2 + 1
    cdelt1 = -1 * celldata / 3600 / oversample
    cdelt2 = celldata / 3600 / oversample    
    modelheader.update('naxis1', nxmod)
    modelheader.update('cdelt1', cdelt1)
    modelheader.update('crpix1', crpix1)
    modelheader.update('crval1', ra_centroid)
    modelheader.update('ctype1', 'RA---SIN')
    modelheader.update('naxis2', nymod)
    modelheader.update('cdelt2', cdelt2)
    modelheader.update('crpix2', crpix2)
    modelheader.update('crval2', dec_centroid)
    modelheader.update('ctype2', 'DEC--SIN')

    g_image, g_lensimage, e_image, e_lensimage, amp_tot, amp_mask = \
            lensutil.sbmap(x, y, nlens, nsource, parameters, model_types)

    #--------------------------------------------------------------------------
    # Compute magnification
    numer = g_lensimage.sum()
    denom = g_image.sum()
    mu_fluxratio = numer / denom
    print 'F_out / F_in = ', mu_fluxratio
    #--------------------------------------------------------------------------

    # write the lensed and unlensed surface brightness maps to disk
    outfits = 'sbmap_image_Region' + sri + '.fits'
    outfits_miriad = 'sbmap_image' + sri
    outfits_unlens = 'sbmap_source' + sri + '.fits'
    os.system('rm -rf ' + outfits)
    os.system('rm -rf ' + outfits_unlens)
    fits.writeto(outfits, g_lensimage, modelheader)
    fits.writeto(outfits_unlens, g_image, modelheader)

    # make the miriad-format model image
    os.system('rm -rf ' + outfits_miriad)
    command = 'fits in=' + outfits + ' op=xyin out=' + outfits_miriad
    os.system(command + ' > dump')


    #--------------------------------------------------------------------------
    # read in Python visibilities
        
    fitsfiles = config.FitsFiles
    #print fitsfiles
    nfiles = len(fitsfiles)

    # read in the observed visibilities
    for file in fitsfiles:

        # read in the observed visibilities
        obsdata = fits.open(file)

        # get the real and imaginary components
        data_real, data_imag, data_wgt = uvutil.visload(obsdata)
        
        #----------------------------------------------------------------------
        # Python version of UVMODEL
        # "Observe" the lensed emission with the SMA and write to a new file
        #----------------------------------------------------------------------
        # Python replace:
        nameindx = file.find('uvfits')
        name = file[0:nameindx-1]

        py_replace = uvmodel.replace(outfits, file)
        #py_real, py_imag, py_wgt = uvmodel.components(py_replace)
        #dreal = py_real - data_real
        #dimag = py_imag - data_imag
        #py_chi2 = py_wgt * (dreal ** 2 + dimag ** 2)
        #print name, 'py-replace ', py_chi2.sum(), py_wgt.sum()

        modelvisfile = name + '_pyreplace.uvfits'
        os.system('rm -rf ' + modelvisfile)
        py_replace.writeto(modelvisfile)

        # and convert to miriad format for imaging
        miriadmodelvisloc = name + '_pyreplace.miriad'
        os.system('rm -rf ' + miriadmodelvisloc)
        command = 'fits in=' + modelvisfile + ' op=uvin out=' + \
            miriadmodelvisloc
        os.system(command + ' > dump')

        # insert fake system temperatures and jy/k values for PdBI data
        if obsdata[0].header['TELESCOP'] == 'PdBI':
            command = 'puthd in=' + miriadmodelvisloc + '/systemp value=40.'
            os.system(command + ' > dump')
            command = 'puthd in=' + miriadmodelvisloc + '/jyperk value=10.'
            os.system(command + ' > dump')
        
        # Python subtract:
        py_subtract = uvmodel.subtract(outfits, file)
        #py_real, py_imag, py_wgt = uvmodel.components(py_subtract)
        #py_chi2 = py_wgt * (py_real ** 2 + py_imag ** 2)
        #print name, 'py-subtract ', py_chi2.sum(), py_wgt.sum()

        modelvisfile = name + '_pysubtract.uvfits'
        os.system('rm -rf ' + modelvisfile)
        py_subtract.writeto(modelvisfile)

        # and convert to miriad format for imaging
        miriadmodelvisloc = name + '_pysubtract.miriad'
        os.system('rm -rf ' + miriadmodelvisloc)
        command = 'fits in=' + modelvisfile + ' op=uvin out=' + \
            miriadmodelvisloc
        os.system(command + ' > dump')

        # insert fake system temperatures and jy/k values for PdBI data
        if obsdata[0].header['TELESCOP'] == 'PdBI':
            command = 'puthd in=' + miriadmodelvisloc + '/systemp value=40.'
            os.system(command + ' > dump')
            command = 'puthd in=' + miriadmodelvisloc + '/jyperk value=10.'
            os.system(command + ' > dump')

    #--------------------------------------------------------------------------
    # Invert the simulated SMA visibilities and deconvolve 

    # the simulated model visibilities
    command = 'csh image_replace.csh'
    os.system(command + ' > dump')

    # the simulated residual visibilities
    command = 'csh image_subtract.csh'
    os.system(command + ' > dump')

    # read in the cleaned images
    simimloc = objectname + '_pyreplace.cm1.fits'
    pythonreplaceim = fits.getdata(simimloc)
    pythonreplaceim = pythonreplaceim[0,0,:,:].copy()

    simimloc = objectname + '_pysubtract.cm1.fits'
    pythonsubtractim = fits.getdata(simimloc)
    pythonsubtractim = pythonsubtractim[0,0,:,:].copy()
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Step 5: Make cutout images and write as fits files

    # default distance from plot axis is 1.05
    beamshift = 1.05

    # get the image centroid in model pixel coordinates
    datawcs = wcs.WCS(headim, naxis=2)
    pix = datawcs.wcs_world2pix(ra_centroid, dec_centroid, 1)
    x0 = pix[0]
    y0 = pix[1]
    modrady = pixextent# nymod / 2.
    modradx = pixextent# nxmod / 2.

    # make data cutout
    totdx1 = x0 - modradx
    totdx2 = x0 + modradx
    totdy1 = y0 - modrady
    totdy2 = y0 + modrady
    imcut = im[totdy1:totdy2,totdx1:totdx2]

    # make sure the image cutouts are square
    #if imcut[0, :].size == nx_contour - 1:
    #    totdx2 = x0 + modradx + 1

    # make cleaned model cutout
    pythoncleancut = pythonreplaceim[totdy1:totdy2,totdx1:totdx2]

    # make cleaned residual cutout
    pythonresidcut = pythonsubtractim[totdy1:totdy2,totdx1:totdx2]

    # make magnification map cutout
    #dmucut = dmu[totdy1:totdy2,totdx1:totdx2]

    fig = plt.figure(figsize=(6.0, 3.0))
    ax = fig.add_subplot(1, 2, 1)

    # overplot position of the background source(s) and the foreground lens(es)
    #nlens = nlens_models[regioni]
    #nsource = nsource_models[regioni]

    for i in range(nsource):
        i6 = i * 6
        xxx = parameters[i6 + 2 + nparlens]
        yyy = parameters[i6 + 3 + nparlens]
        source_pa = 90 - parameters[i6 + 5 + nparlens]
        model_type = model_types[i]
        if model_type == 'gaussian':
            norm = 2.35
        if model_type == 'cylinder':
            norm = numpy.sqrt(2)
        meansize = norm * parameters[i6 + 1 + nparlens]
        source_bmaj = meansize / numpy.sqrt(parameters[i6 + 4 + nparlens])
        source_bmin = meansize * numpy.sqrt(parameters[i6 + 4 + nparlens])
        e = Ellipse((xxx, yyy), source_bmaj, source_bmin, \
                angle=source_pa, ec='white', lw=0.5, fc='magenta', \
                zorder=2, fill=True, alpha=0.5)
        ax.add_artist(e)
    for i in range(nlens):
        i5 = i * 5
        xxx = numpy.array([parameters[i5 + 1]])
        yyy = numpy.array([parameters[i5 + 2]])
        plt.plot(xxx, yyy, 'o', ms=5., mfc='black', mec='white', mew=0.5, \
            label='Lens Position', zorder=20)
        lens_pa = 90 - parameters[i5 + 4]
        meansize = 2 * parameters[i5]
        lens_bmaj = meansize / numpy.sqrt(parameters[i5 + 3])
        lens_bmin = meansize * numpy.sqrt(parameters[i5 + 3])
        elens = Ellipse((xxx, yyy), lens_bmaj, lens_bmin, \
                angle=lens_pa, ec='orange', lw=1.0, \
                zorder=20, fill=False)
        ax.add_artist(elens)
        #plt.plot(-xxx, yyy, '+', ms=3., color='white')


    ndof = (2*modradx+1) * (2*modrady+1)

    plt.subplot(1, 2, 1)

    cellplus = celldata*(2*modradx+1.1)/(2*modradx)
    modx = ( numpy.arange(2*modradx) - modradx ) * (-cellplus) - celldata/2.
    mody = ( numpy.arange(2*modradx) - modradx ) * cellplus + celldata/2.
    extent = [modx[0], modx[-1], mody[0], mody[-1] ]
    plt.imshow(pythoncleancut, cmap='gray_r', interpolation='nearest', extent=extent, \
            origin='lower')

    plevs = rms * 2**(numpy.arange(10) + 1)
    nlevs = sorted(-rms * 2**(numpy.arange(4) + 1))
    pcline = 'solid'
    ncline = 'dashed'
    nx_contour = imcut[0, :].size
    ny_contour = imcut[:, 0].size
    cmodx = ( numpy.arange(nx_contour) - modradx ) * (-cellplus) - celldata/2.
    cmody = ( numpy.arange(ny_contour) - modradx ) * cellplus + celldata/2.
    plt.contour(cmodx, cmody, imcut, colors='red', levels=plevs, linestyles=pcline, \
        linewidths=1.5)
    plt.contour(cmodx, cmody, imcut, colors='red', levels=nlevs, linestyles=ncline, \
        linewidths=1.5)

    # plot the critical curve
    #plt.contour(cmodx, cmody, dmu, colors='orange', levels=[100])

    plt.minorticks_on()
    plt.tick_params(width=1.5, which='both')
    plt.tick_params(length=2, which='minor')
    plt.tick_params(length=4, which='major')
    plt.xlabel(r'$\Delta$RA (arcsec)', fontsize='x-large')
    plt.ylabel(r'$\Delta$Dec (arcsec)', fontsize='x-large')

    axisrange = [numpy.float(xhi),numpy.float(xlo),numpy.float(ylo),numpy.float(yhi)]
    plt.axis(axisrange)

    normx = axisrange[0] - axisrange[1]
    normy = axisrange[3] - axisrange[2]
    bparad = bpa / 180 * math.pi
    beamx = numpy.abs(numpy.sin(bparad) * bmaj) + numpy.abs(numpy.cos(bparad) * bmin)
    beamy = numpy.abs(numpy.cos(bparad) * bmaj) + numpy.abs(numpy.sin(bparad) * bmin)
    dx = numpy.float(xhi) - numpy.float(xlo)
    dy = numpy.float(yhi) - numpy.float(ylo)
    bufferx = 0.03 * dx / 6.0
    buffery = 0.03 * dx / 6.0
    xpos = 1 - beamx/dx/2 - bufferx
    ypos = beamy/dy/2 + buffery
    #beamx = bmaj * numpy.abs(numpy.cos(bparad))
    #beamy = bmaj * numpy.abs(numpy.sin(bparad))
    xpos = 0.95 * axisrange[1] + 0.95 * beamx / 2.
    ypos = 0.95 * axisrange[2] + 0.95 * beamy / 2.

    e = Ellipse((xpos,ypos), bmaj, bmin, angle=90 - bpa, ec='black', \
        hatch='///', lw=1.5, fc='None', zorder=10, fill=True)
    ax.add_artist(e)
    #plt.text(0.08, 0.85, objname, transform=ax.transAxes, fontsize='x-large')

    objtitle = objectname + '.' + str(regioni + 1)
    plt.suptitle(objtitle, fontsize='x-large', y=0.99)

    # plot the difference between the cleaned image and the sma data
    plt.subplot(1, 2, 2)
    ax1 = fig.add_subplot(1, 2, 2)

    modx = ( numpy.arange(2*modradx) - modradx ) * (-cellplus) - celldata/2.
    mody = ( numpy.arange(2*modrady) - modrady ) * cellplus + celldata/2.
    nrms = -5*rms
    prms = 5*rms
    extent = [modx[0], modx[-1], mody[0], mody[-1] ]
    plt.imshow(pythonresidcut, cmap='gray_r', interpolation='nearest', extent=extent, \
            origin='lower', vmax=3*rms, vmin = -3*rms)

    plevs = rms * 2**(numpy.arange(4) + 1)
    nlevs = sorted(-rms * 2**(numpy.arange(4) + 1))
    pcline = 'solid'
    ncline = 'dashed'
    cmodx = ( numpy.arange(nx_contour) - modradx ) * (-cellplus) - celldata/2.
    cmody = ( numpy.arange(ny_contour) - modradx ) * cellplus + celldata/2.
    plt.contour(cmodx, cmody, pythonresidcut, colors='white', levels=plevs,
        linewidths=1.5, cline=pcline)
    plt.contour(cmodx, cmody, pythonresidcut, colors='black', levels=nlevs,
        linewidths=1.5, cline=ncline)

    for i in range(nsource):
        i6 = i * 6
        xxx = parameters[i6 + 2 + nparlens]
        yyy = parameters[i6 + 3 + nparlens]
        source_pa = 90 - parameters[i6 + 5 + nparlens]
        model_type = model_types[i]
        if model_type == 'gaussian':
            norm = 2.35
        if model_type == 'cylinder':
            norm = numpy.sqrt(2)
        meansize = norm * parameters[i6 + 1 + nparlens]
        source_bmaj = meansize / numpy.sqrt(parameters[i6 + 4 + nparlens])
        source_bmin = meansize * numpy.sqrt(parameters[i6 + 4 + nparlens])
        e = Ellipse((xxx, yyy), source_bmaj, source_bmin, \
                angle=source_pa, ec='white', lw=0.5, fc='magenta', \
                zorder=2, fill=True, alpha=0.5)
        ax1.add_artist(e)
    for i in range(nlens):
        i5 = i * 5
        xxx = numpy.array([parameters[i5 + 1]])
        yyy = numpy.array([parameters[i5 + 2]])
        plt.plot(xxx, yyy, 'o', ms=5., mfc='black', mec='white', mew=0.5, \
            label='Lens Position', zorder=20)
        lens_pa = 90 - parameters[i5 + 4]
        meansize = 2 * parameters[i5]
        lens_bmaj = meansize / numpy.sqrt(parameters[i5 + 3])
        lens_bmin = meansize * numpy.sqrt(parameters[i5 + 3])
        elens = Ellipse((xxx, yyy), lens_bmaj, lens_bmin, \
                angle=lens_pa, ec='orange', lw=1.0, \
                zorder=20, fill=False)
        ax1.add_artist(elens)
        #plt.plot(-xxx, yyy, '+', ms=3., color='white')

    plt.minorticks_on()
    plt.tick_params(width=1.5, which='both')
    plt.tick_params(length=2, which='minor')
    plt.tick_params(length=4, which='major')
    plt.xlabel(r'$\Delta$RA (arcsec)', fontsize='x-large')
    plt.ylabel(r'$\Delta$Dec (arcsec)', fontsize='x-large')

    plt.axis([numpy.float(xhi), numpy.float(xlo), numpy.float(ylo),
        numpy.float(yhi)])

    e2 = Ellipse((xpos,ypos), bmaj, bmin, angle=90 - bpa, ec='black', \
        hatch='///', lw=1.5, fc='None', zorder=10, fill=True)
    ax1.add_artist(e2)

    #plt.tight_layout(pad=1.5)
    plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.10, wspace=0.35)
    #ax1.set_aspect('equal')

    # rename files to remove '.' from the name
    rename = objectname.replace('.', '_')
    #os.chdir(cwd)
    #savefig('sbmap_purepython/' + rename + '_sbmap')
    savefig(rename + '.sbmap.Region' + sri + '.pdf')

cmd = 'rm -rf sbmap* *pyreplace* *pysubtract* dump'
os.system(cmd)
