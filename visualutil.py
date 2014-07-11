
def plotPDF(fitresults, tag, limits='', Ngood=5000, axes='auto'):

    """

    Plot the PDF of each parameter of the model.

    """

    import numpy
    import matplotlib.pyplot as plt
    from pylab import savefig
    from matplotlib import rc
    import modifypdf


    # plotting parameters
    rc('font',**{'family':'sans-serif', 'sans-serif':['Arial Narrow'], 
        'size':'12'})

    # grab the last Ngood fits
    fitresults = fitresults[-Ngood:]
    lnprobstring = "prior to pruning <Ln Prob>: {:f}"
    #print(lnprobstring.format(fitresults['lnprob'].mean()))

    # identify the good fits
    fitresultsgood = modifypdf.prune(fitresults)

    # determine dimensions of PDF plots
    nparams = len(fitresultsgood[0])
    ncol = 4
    nrow = nparams / ncol + 1
    j = 1

    fig = plt.figure(figsize=(12.0, 1.5 * nrow))

    # set up the plotting window
    plt.subplots_adjust(left=0.08, bottom=0.1, right=0.95, top=0.95,
        wspace=0.4, hspace=0.65)

    pnames = fitresultsgood.names

    counter = 0
    for pname in pnames:

        # Position of primary lens
        frg = fitresultsgood[pname]
        rmsval = numpy.std(frg)
        if rmsval > 1e-6:
            avgval = numpy.mean(frg)
            print(pname + ' = ', avgval, ' +/- ', rmsval)
            totalwidth = frg.max() - frg.min()
            nbins = totalwidth / rmsval * 5
            ax = plt.subplot(nrow, ncol, j)
            j += 1
            plt.hist(frg, nbins, edgecolor='blue')
            plt.ylabel('N')
            plt.xlabel(pname)
            if axes == 'auto':
                start, end = ax.get_xlim()
                nticks = 5
                stepsize = (end - start) / nticks
                ax.xaxis.set_ticks(numpy.arange(start, end + 0.99*stepsize, 
                    stepsize))
            elif axes == 'initial':
                oldaxis = plt.axis()
                if pname[0:6] == 'lnprob':
                    xmin = frg.min()
                    xmax = frg.max()
                elif pname[0:2] == 'mu':
                    xmin = 0
                    xmax = 30
                else:
                    p_l = limits[0]
                    p_u = limits[1]
                    xmin = p_l[counter]
                    xmax = p_u[counter]
                    counter += 1
                ymin = oldaxis[2]
                ymax = oldaxis[3]
                plt.axis([xmin, xmax, ymin, ymax])


    savefile = tag + 'PDFs.png'
    savefig(savefile)

def makeSBmap(paramData, parameters, regioni):

    """

    Make a surface brightness map of the lensed image for a given set of model
    parameters.

    """

    import lensutil
    from astropy.io import fits
    import os


    nlens = paramData['nlens_regions'][regioni]
    nsource = paramData['nsource_regions'][regioni]
    x = paramData['x'][regioni]
    y = paramData['y'][regioni]
    modelheader = paramData['modelheader'][regioni]
    model_types = paramData['model_types'][regioni]

    SBmap, LensedSBmap, Aperture, LensedAperture, mu_tot, mu_mask = \
            lensutil.sbmap(x, y, nlens, nsource, parameters, model_types)

    sri = str(regioni)
    LensedSBmapLoc = 'LensedSBmap_Region' + sri + '.fits'
    SBmapLoc = 'SBmap_Region' + sri + '.fits'
    cmd = 'rm -rf ' + LensedSBmapLoc + ' ' + SBmapLoc
    os.system(cmd)

    fits.writeto(LensedSBmapLoc, LensedSBmap, modelheader)
    fits.writeto(SBmapLoc, SBmap, modelheader)

    return

def makeVis(config, regioni):

    """

    Make simulated visibilities given a model image and observed visibilities.
    Writes the visibilities to uvfits files.

    """

    import uvmodel
    import os


    visfile = config['UVData']
    
    #----------------------------------------------------------------------
    # Python version of UVMODEL
    # "Observe" the lensed emission with the SMA and write to a new file
    #----------------------------------------------------------------------
    # Python version of UVMODEL's "replace" subroutine:
    nameindx = visfile.find('uvfits')
    name = visfile[0:nameindx-1]
    print(name)
    visfile = name + '.ms'

    sri = str(regioni)
    SBmapLoc = 'LensedSBmap_Region' + sri + '.fits'
    modelvisfile = name + '.Region' + sri + '_model.ms'
    os.system('rm -rf ' + modelvisfile)
    uvmodel.replace(SBmapLoc, visfile, modelvisfile)
    
    # Python version of UVMODEL's "subtract" subroutine:
    modelvisfile = name + '.Region' + sri + '_residual.ms'
    os.system('rm -rf ' + modelvisfile)
    uvmodel.subtract(SBmapLoc, visfile, modelvisfile)

def makeImage(config, objectname, regioni, interactive=True, miriad=False):

    """

    Make an image of the model and the residual from simulated model
    visibilities.  Requires CASA or miriad.

    """

    import os
    from astropy.io import fits

        
    visfile = config['UVData']
    target = config['ObjectName']
    sri = str(regioni)
    fitsim = config['ImageName']
    fitshead = fits.getheader(fitsim)
    imsize = [fitshead['NAXIS1'], fitshead['NAXIS2']]
    cell = str(fitshead['CDELT2'] * 3600) + 'arcsec'

    # search for an existing mask
    index = fitsim.find('.fits')
    maskname = fitsim[0:index] + '.mask'
    try:
        maskcheck = os.path.exists(maskname)
    except:
        maskcheck = False

    if maskcheck:
        mask = maskname
    else:
        mask = ''

    # invert and clean the simulated model visibilities
    imloc = target + '_Region' + sri + '_model'

    if miriad:
        # use miriad for imaging
        index = visfile.find('.uvfits')
        name = visfile[0:index]
        miriadmodelvisloc = name + '.Region' + sri + '_model.miriad'
        miriadin = miriadmodelvisloc + ' ' + imloc + ' ' + imsize + ' ' + cell
        command = 'csh image.csh ' + miriadin
        os.system(command + ' > dump')

        # the simulated residual visibilities
        imloc = target + '_Region' + sri + '_residual'
        miriadresidvisloc = name + '.Region' + sri + '_residual.miriad'
        miriadin = miriadresidvisloc + ' ' + imloc + ' ' + imsize + ' ' + cell
        command = 'csh image.csh ' + miriadin
        os.system(command + ' > dump')
    else:
        # use CASA for imaging
        from clean import clean
        from casa import exportfits
        index = visfile.find('.ms')
        name = visfile[0:index]
        miriadmodelvisloc = name + '.Region' + sri + '_model.ms'
        os.system('rm -rf ' + imloc + '*')
        clean(vis=miriadmodelvisloc, imagename=imloc, mode='mfs', niter=10000,
            threshold='0.2mJy', interactive=interactive, mask=mask, 
            imsize=imsize, cell=cell, weighting='briggs', robust=0.5)

        # export the cleaned image to a fits file
        os.system('rm -rf ' + imloc + '.fits')
        exportfits(imagename=imloc + '.image', fitsimage=imloc + '.fits')

        # invert and clean the residual visibilities
        imloc = target + '_Region' + sri + '_residual'
        miriadmodelvisloc = name + '.Region' + sri + '_residual.ms'
        os.system('rm -rf ' + imloc + '*')
        clean(vis=miriadmodelvisloc, imagename=imloc, mode='mfs', niter=10000,
            threshold='0.2mJy', interactive=interactive, mask=mask, 
            imsize=imsize, cell=cell, weighting='briggs', robust=0.5)

        # export the cleaned image to a fits file
        os.system('rm -rf ' + imloc + '.fits')
        exportfits(imagename=imloc + '.image', fitsimage=imloc + '.fits')
    
    return

def plotImage(model, data, config, parameters, regioni, modeltype, tag=''):

    """

    Make a surface brightness map of a given model image.  Overlay with red
    contours a surface brightness map of the data image to which the model was
    fit. 

    """

    import numpy
    from astropy import wcs
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    from pylab import savefig
    import setuputil
    import re

    # set font properties
    font = {'family' : 'Arial Narrow',
            'weight' : 'bold',
            'size'   : 10}
    matplotlib.rc('font', **font)
    matplotlib.rcParams['axes.linewidth'] = 1.5

    fig = plt.figure(figsize=(3.0, 3.0))
    ax = fig.add_subplot(1, 1, 1)
    plt.subplots_adjust(left=0.08, right=0.97, top=0.97, 
            bottom=0.08, wspace=0.35)

    configkeys = config.keys()
    configkeystring = " ".join(configkeys)
    regionlist = re.findall('Region.', configkeystring)
    region = regionlist[regioni]
    cr = config[region]

    ra_centroid = cr['RACentroid']
    dec_centroid = cr['DecCentroid']
    radialextent = cr['RadialExtent']

    # count the number of lenses
    configkeys = cr.keys()
    configkeystring = " ".join(configkeys)
    lenslist = re.findall('Lens.', configkeystring)
    nlens = len(lenslist)

    # count the number of sources
    sourcelist = re.findall('Source.', configkeystring)
    nsource = len(sourcelist)

    nparlens = 5 * nlens

    # get the image centroid in model pixel coordinates
    headim = data[0].header
    headmod = model[0].header
    im = data[0].data
    im = im[0, 0, :, :]

    # good region is where mask is zero
    mask = setuputil.makeMask(config)
    goodregion = mask == 0

    # compute sigma image from cutout of SMA flux image
    rms = im[goodregion].std()

    # Obtain measurements of beamsize and image min/max
    bmaj = headim['BMAJ'] * 3600
    bmin = headim['BMIN'] * 3600
    bpa = headim['BPA']
    cdelt1 = headim['CDELT1'] * 3600
    cdelt2 = headim['CDELT2'] * 3600
    cell = numpy.sqrt( abs(cdelt1) * abs(cdelt2) )

    im_model = model[0].data
    if im_model.ndim == 4:
        im_model = im_model[0, 0, :, :]
    #nx_model = im_model[0, :].size
    pixextent = radialextent / cell
    datawcs = wcs.WCS(headim, naxis=2)
    pix = datawcs.wcs_world2pix(ra_centroid, dec_centroid, 1)
    x0 = numpy.round(pix[0])
    y0 = numpy.round(pix[1])
    imrady = radialextent / cell# nymod / 2.
    imradx = radialextent / cell# nxmod / 2.

    # make data cutout
    totdx1 = x0 - imradx
    totdx2 = x0 + imradx
    totdy1 = y0 - imrady
    totdy2 = y0 + imrady
    datacut = im[totdy1:totdy2,totdx1:totdx2]

    # make cleaned model cutout
    headerkeys = headmod.keys()
    cd1_1 = headerkeys.count('CD1_1')
    cd1_2 = headerkeys.count('CD1_2')
    if cd1_1 == 0:
        cdelt1_model = numpy.abs(headmod['CDELT1'] * 3600)
        cdelt2_model = numpy.abs(headmod['CDELT2'] * 3600)
    else:
        cdelt1_model = numpy.abs(headmod['CD1_1'] * 3600)
        cdelt2_model = numpy.abs(headmod['CD2_2'] * 3600)
        cd11 = headmod['CD1_1']
        if cd1_2 == 0:
            cd12 = 0
            cd21 = 0
        else:
            cd12 = headmod['CD1_2']
            cd21 = headmod['CD2_1']
        cd22 = headmod['CD2_2']
        cdelt1_model = numpy.sqrt(cd11 ** 2 + cd12 ** 2) * 3600
        cdelt2_model = numpy.sqrt(cd21 ** 2 + cd22 ** 2) * 3600
        if cd12 == 0:
            cd12 = cd11 / 1e8
        cdratio = numpy.abs(cd11 / cd12)
        if cdratio < 1:
            cdratio = 1 / cdratio
    cellmod = numpy.sqrt( abs(cdelt1_model) * abs(cdelt2_model) )

    modelwcs = wcs.WCS(headmod, naxis=2)
    pix = modelwcs.wcs_world2pix(ra_centroid, dec_centroid, 1)
    x0 = numpy.round(pix[0])
    y0 = numpy.round(pix[1])
    modrady = radialextent / cellmod
    modradx = radialextent / cellmod
    totdx1 = x0 - modradx
    totdx2 = x0 + modradx
    totdy1 = y0 - modrady
    totdy2 = y0 + modrady
    modelcut = im_model[totdy1:totdy2,totdx1:totdx2]

    for i in range(nsource):
        i6 = i * 6
        xxx = parameters[i6 + 2 + nparlens]
        yyy = parameters[i6 + 3 + nparlens]
        source_pa = 90 - parameters[i6 + 5 + nparlens]
        #model_type = model_types[i]
        #if model_type == 'gaussian':
        norm = 2.35
        #if model_type == 'cylinder':
        #    norm = numpy.sqrt(2)
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

    #cellp = cell * (2 * pixextent + 1.1) / (2 * pixextent)
    xlo = -radialextent
    xhi = radialextent
    ylo = -radialextent
    yhi = radialextent
    ncell = (xhi - xlo) / cell
    modx = -numpy.linspace(xlo, xhi, ncell)
    mody = numpy.linspace(ylo, yhi, ncell)
    #modx = -(numpy.arange(2 * pixextent) - pixextent) * cellp - cell/2.
    #mody = (numpy.arange(2 * pixextent) - pixextent) * cellp + cell/2.
    cornerextent = [modx[0], modx[-1], mody[0], mody[-1] ]
    if modeltype == 'residual':
        grayscalename = 'Residual'
        pcolor = 'white'
        ncolor = 'black'
        vmax = 5 * rms
        vmin = -5 * rms
    elif modeltype == 'model':
        grayscalename = 'Model'
        pcolor = 'red'
        ncolor = 'red'
        vmax = modelcut.max()
        vmin = -3 * rms
    else:
        grayscalename = config['OpticalTag']
        filtindx = grayscalename.find(' ')
        filtname = grayscalename[filtindx + 1:]
        if filtname == 'F110W':
            modelcut = numpy.log10(modelcut - modelcut.min() + 1)
        pcolor = 'red'
        ncolor = 'red'
        vmax = modelcut.max()
        vmin = modelcut.min()
    plt.imshow(modelcut, cmap='gray_r', interpolation='nearest', \
            extent=cornerextent, origin='lower', vmax=vmax, vmin=vmin)

    plevs = 3*rms * 2**(numpy.arange(10))
    nlevs = sorted(-3 * rms * 2**(numpy.arange(4)))
    pcline = 'solid'
    ncline = 'dashed'
    #nx_contour = datacut[0, :].size
    #ny_contour = datacut[:, 0].size
    #cmodx = -(numpy.arange(nx_contour) - pixextent) * cellp - cell/2.
    #cmody = (numpy.arange(ny_contour) - pixextent) * cellp + cell/2.
    cmodx = -numpy.linspace(xlo, xhi, ncell)
    cmody = numpy.linspace(ylo, yhi, ncell)
    plt.contour(cmodx, cmody, datacut, colors=pcolor, levels=plevs, \
            linestyles=pcline, linewidths=1.5)
    plt.contour(cmodx, cmody, datacut, colors=ncolor, levels=nlevs, \
            linestyles=ncline, linewidths=1.5)

    # plot the critical curve
    #plt.contour(cmodx, cmody, dmu, colors='orange', levels=[100])
    #axisrange = plt.axis()
    axisrange = [numpy.float(xhi),numpy.float(xlo),numpy.float(ylo),numpy.float(yhi)]
    plt.axis(axisrange)

    plt.minorticks_on()
    plt.tick_params(width=1.5, which='both')
    plt.tick_params(length=2, which='minor')
    plt.tick_params(length=4, which='major')
    #plt.xlabel(r'$\Delta$RA (arcsec)', fontsize='x-large')
    #plt.ylabel(r'$\Delta$Dec (arcsec)', fontsize='x-large')

    bparad = bpa / 180 * numpy.pi
    beamx = numpy.abs(numpy.sin(bparad) * bmaj) + \
            numpy.abs(numpy.cos(bparad) * bmin)
    beamy = numpy.abs(numpy.cos(bparad) * bmaj) + \
            numpy.abs(numpy.sin(bparad) * bmin)
    beamxhi = 2 * pixextent / cell
    beamxlo = -2 * pixextent / cell
    beamyhi = 2 * pixextent / cell
    beamylo = -2 * pixextent / cell
    beamdx = numpy.float(beamxhi) - numpy.float(beamxlo)
    beamdy = numpy.float(beamyhi) - numpy.float(beamylo)
    bufferx = 0.03 * beamdx / 6.0
    buffery = 0.03 * beamdx / 6.0
    xpos = 1 - beamx/beamdx/2 - bufferx
    ypos = beamy/beamdy/2 + buffery
    #beamx = bmaj * numpy.abs(numpy.cos(bparad))
    #beamy = bmaj * numpy.abs(numpy.sin(bparad))

    xpos = 0.95 * axisrange[1] + 0.95 * beamx / 2.
    ypos = 0.95 * axisrange[2] + 0.95 * beamy / 2.

    e = Ellipse((xpos,ypos), bmaj, bmin, angle=90 - bpa, ec='black', \
        hatch='//////', lw=1.0, fc='None', zorder=10, fill=True)
    ax.add_artist(e)

    plt.text(0.92, 0.88, grayscalename, transform=ax.transAxes,
            fontsize='xx-large', ha='right')
    objname = config['ObjectName']
    plt.text(0.08, 0.88, objname, transform=ax.transAxes,
            fontsize='xx-large')
    sri = str(regioni)
    bigtag = '.' + modeltype + '.' + tag

    savefig('LensedSBmap.Region' + sri + bigtag + '.pdf')
    #plt.clf()

def removeTempFiles():

    """

    Remove files created along the way by visualutil routines.

    """

    import os


    cmd = 'rm -rf *SBmap*fits *_model* *_residual* dump'
    os.system(cmd)    

def plotFit(config, paramData, parameters, regioni, tag='', cleanup=True,
        showOptical=False, interactive=True):

    """

    Plot a particular model fit.

    """

    from astropy.io import fits


    # make the lensed image
    makeSBmap(paramData, parameters, regioni)

    # make the simulated visibilities
    makeVis(config, regioni)

    # image the simulated visibilities
    objectname = config['ObjectName']
    makeImage(config, objectname, regioni, interactive=interactive)

    # read in the images of the simulated visibilities
    sri = str(regioni)
    simimloc = objectname + '_Region' + sri + '_model.fits'
    model = fits.open(simimloc)

    simimloc = objectname + '_Region' + sri + '_residual.fits'
    residual = fits.open(simimloc)

    # read in the data
    data = fits.open(config['ImageName'])

    if showOptical:
        # read in the data
        optical = fits.open(config['OpticalImage'])

        # plot the images
        plotImage(optical, data, config, parameters, regioni, 'optical', 
                tag=tag)

    # plot the images
    plotImage(model, data, config, parameters, regioni, 'model', tag=tag)

    # plot the residual
    plotImage(residual, residual, config, parameters, regioni, 'residual', 
            tag=tag)

    # remove the intermediate files
    if cleanup:
        removeTempFiles()

def preProcess(config, paramData, fitresult, tag='', cleanup=True,
        showOptical=False, interactive=True):

    """

    Cycle through each region and run plotFit, selecting parameters
    appropriately.

    """

    import setuputil
    import numpy
    import re


    # Loop over each region
    nlensedsource = paramData['nlensedsource']
    nlensedregions = paramData['nlensedregions']
    npar_previous = 0

    configkeys = config.keys()
    configkeystring = " ".join(configkeys)
    regionlist = re.findall('Region.', configkeystring)
    for regioni, region in enumerate(regionlist):
        cr = config[region]

        nmu = 0#2 * (numpy.array(nlensedsource).sum() + nlensedregions)
        if nmu > 0:
            allparameters0 = list(fitresult)[1:-nmu]
        else:
            allparameters0 = list(fitresult)[1:]

        # search poff_models for parameters fixed relative to other parameters
        fixindx = paramData['fixindx']
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
                par0 = fitresult[fixindx[subindx] + 1]
            parameters_offset[ifixed] = fitresult[subindx + 1] + par0

        allparameters = allparameters0 + parameters_offset

        # count the number of lenses
        configkeys = cr.keys()
        configkeystring = " ".join(configkeys)
        lenslist = re.findall('Lens.', configkeystring)
        nlens = len(lenslist)

        # count the number of sources
        sourcelist = re.findall('Source.', configkeystring)
        nsource = len(sourcelist)

        nparlens = 5 * nlens
        nparsource = 6 * nsource
        npar = nparlens + nparsource + npar_previous
        parameters = allparameters[npar_previous:npar]
        npar_previous = npar
        plotFit(config, paramData, parameters, regioni, tag=tag,
                cleanup=cleanup, showOptical=showOptical,
                interactive=interactive)
