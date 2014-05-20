import numpy
from astropy.io import fits


def loadParams(config):
    headim = fits.getheader(config.ImageName)

    # get resolution in ALMA image
    celldata = numpy.abs(headim['CDELT1'] * 3600)

    #--------------------------------------------------------------------------
    # Define the number of walkers
    nwalkers = config.Nwalkers

    # Determine method of computing lnlike
    lnlikemethod = config.lnLike

    # determine the number of regions for which we need surface brightness maps
    regionIDs = config.RegionID
    nregions = len(regionIDs)

    # instantiate lists that must be carried through to lnprob function
    x = []
    y = []
    modelheader = []
    nlens_regions = []
    nsource_regions = []
    p_u = []
    p_l = []
    prior_shape = []
    poff = []
    pname = []
    pzero = []
    model_types = []
    nparams_total = 0
    nlensedsource = 0
    nlensedregions = 0
    for i in range(nregions):
        ri = str(i)
        ra_centroid = config.RACentroid[i]
        dec_centroid = config.DecCentroid[i]
        extent = config.RadialExtent[i]
        oversample = config.Oversample[i]
        nlens = config.Nlens[i]
        nsource = config.Nsource[i]

        # Append the number of lenses and sources for this region
        nlens_regions.append(nlens)
        nsource_regions.append(nsource)

        # define number of pixels in lensed surface brightness map
        dx = 2 * extent
        nxmod = oversample * int(round(dx / celldata))
        dy = 2 * extent
        nymod = oversample * int(round(dy / celldata))

        # make x and y coordinate images for lens model
        onex = numpy.ones(nxmod)
        oney = numpy.ones(nymod)
        linspacex = numpy.linspace(0, 1, nxmod)
        linspacey = numpy.linspace(0, 1, nymod)
        x.append(dx * numpy.outer(oney, linspacex) - extent)
        y.append(dy * numpy.outer(linspacey, onex) - extent)

        # Provide world-coordinate system transformation data in the header of
        # the lensed surface brightness map
        headmod = headim.copy()
        crpix1 = nxmod / 2 + 1
        crpix2 = nymod / 2 + 1
        cdelt1 = -1 * celldata / 3600 / oversample
        cdelt2 = celldata / 3600 / oversample
        headmod['naxis1'] = nxmod
        headmod['cdelt1'] = cdelt1
        headmod['crpix1'] = crpix1
        headmod['crval1'] = ra_centroid
        headmod['ctype1'] = 'RA---SIN'
        headmod['naxis2'] = nymod
        headmod['cdelt2'] = cdelt2
        headmod['crpix2'] = crpix2
        headmod['crval2'] = dec_centroid
        headmod['ctype2'] = 'DEC--SIN'
        modelheader.append(headmod)

        # the parameter initialization vectors
        p1 = []
        p2 = []

        for ilens in range(nlens):

            li = str(ilens)

            # constraints on the lenses
            lensparams = ['EinsteinRadius', 
                    'DeltaRA', 
                    'DeltaDec', 
                    'AxialRatio',
                    'PositionAngle']
            tag = '_Lens' + li + '_Region' + ri
            for lensparam in lensparams:
                fullparname = 'Prior_' + lensparam + tag
                values = getattr(config, fullparname)
                prior_shape.append(values.pop()) 
                #values = getattr(config, fullparname)
                poff.append(values.pop()) 
                values = numpy.array(values).astype(float)
                if values.size < 2:
                    import pdb; pdb.set_trace()
                p_u.append(values[1]) 
                p_l.append(values[0]) 
                pname.append(lensparam + tag)
                fullparname = 'Init_' + lensparam + tag
                values = getattr(config, fullparname)
                p2.append(values[1]) 
                p1.append(values[0]) 

        model_types_source = []
        if nlens > 0:
            nlensedsource += nsource
            nlensedregions += 1
        for isource in range(nsource):

            si = str(isource)

            sourceparams = ['IntrinsicFlux', 
                    'Size', 
                    'DeltaRA', 
                    'DeltaDec',
                    'AxialRatio', 
                    'PositionAngle']
            tag = '_Source' + si + '_Region' + ri
            for sourceparam in sourceparams:
                fullparname = 'Prior_' + sourceparam + tag
                values = getattr(config, fullparname)
                prior_shape.append(values.pop()) 
                #values = getattr(config, fullparname)
                poff.append(values.pop()) 
                values = numpy.array(values).astype(float)
                p_u.append(values[1]) 
                p_l.append(values[0]) 
                pname.append(sourceparam + tag)
                fullparname = 'Init_' + sourceparam + tag
                values = getattr(config, fullparname)
                p2.append(values[1]) 
                p1.append(values[0]) 

            # get the model type
            fullparname = 'ModelMorphology' + tag
            model_types_source.append(getattr(config, fullparname))

        # append the set of model types for this region
        model_types.append(model_types_source)

        # determine the number of free parameters in the model
        nparams = len(p1)

        # add that number to the total number of free parameters considering
        # all regions so far
        nparams_total += nparams

        # Otherwise, choose an initial set of positions for the walkers.
        pzero_model = numpy.zeros((nwalkers, nparams))
        for j in range(nparams):
            #if p3[j] == 'uniform':
            pzero_model[:, j] = numpy.random.uniform(p1[j], p2[j], nwalkers)
            #if p3[j] == 'normal':
            #    pzero_model[:,j] = (numpy.random.normal(loc=p1[j], 
            #    scale=p2[j], size=nwalkers))
            #if p4[j] == 'pos':
            #    pzero[:, j] = numpy.abs(pzero[:, j])
        if pzero == []:
            pzero = pzero_model
        else:
            pzero = numpy.append(pzero, pzero_model, axis=1)

    paramSetup = {'x': x, 
            'y': y, 
            'modelheader': modelheader,
            'nlens_regions': nlens_regions, 
            'nsource_regions': nsource_regions,
            'nlensedsource': nlensedsource,
            'nlensedregions': nlensedregions,
            'p_u': numpy.array(p_u), 
            'p_l': numpy.array(p_l), 
            'prior_shape': numpy.array(prior_shape),
            'poff': poff, 
            'pname': pname, 
            'pzero': pzero, 
            'model_types': model_types, 
            'nwalkers': nwalkers, 
            'nparams': nparams_total, 
            'celldata': celldata,
            'lnlikemethod': lnlikemethod,
            'nregions': nregions}
    return paramSetup

def fixParams(paramSetup):
    """
    Determine the indices for fixed parameters.
    """
    nparams = paramSetup['nparams']
    poff = paramSetup['poff']
    pname = paramSetup['pname']
    fixindx = numpy.zeros(nparams) - 1
    for ifix in range(nparams):
        if pname.count(poff[ifix]) > 0:
            fixindx[ifix] = pname.index(poff[ifix])
    return fixindx

def getCell(headim):
    celldata = numpy.abs(headim['CDELT1'] * 3600)
    return celldata

def makeMask(config):

    from astropy import wcs

    imloc = config.ImageName
    headim = fits.getheader(imloc)
    im = fits.getdata(imloc)
    im = im[0, 0, :, :]
    celldata = getCell(headim)
    datawcs = wcs.WCS(headim, naxis=2)
    nx = im[0,:].size
    ny = im[:,0].size

    # compute rms within central 3/4 of the image
    mask = im.copy()
    mask[:] = 1
    yr0 = ny / 4
    yr1 = 3 * ny / 4
    xr0 = nx / 4
    xr1 = 3 * nx / 4
    mask[yr0:yr1, xr0:xr1] = 0
    nregions = len(config.RACentroid)
    for regioni in range(nregions):
        ra_centroid = config.RACentroid[regioni]
        dec_centroid = config.DecCentroid[regioni]
        extent = config.RadialExtent[regioni]

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
    return mask 
