#
# lensdemo_funcs.py
#
# Function module for strong lensing demos
#
# Intended for use with lensdemo_script.py
#
# Copyright 2009 by Adam S. Bolton
# Creative Commons Attribution-Noncommercial-ShareAlike 3.0 license applies:
# http://creativecommons.org/licenses/by-nc-sa/3.0/
# All redistributions, modified or otherwise, must include this
# original copyright notice, licensing statement, and disclaimer.
# DISCLAIMER: ABSOLUTELY NO WARRANTY EXPRESS OR IMPLIED.
# AUTHOR ASSUMES NO LIABILITY IN CONNECTION WITH THIS COMPUTER CODE.
#

import numpy as N

def xy_rotate(x, y, xcen, ycen, phi):
    """
    NAME: xy_rotate

    PURPOSE: Transform input (x, y) coordiantes into the frame of a new
             (x, y) coordinate system that has its origin at the point
             (xcen, ycen) in the old system, and whose x-axis is rotated
             c.c.w. by phi degrees with respect to the original x axis.

    USAGE: (xnew,ynew) = xy_rotate(x, y, xcen, ycen, phi)

    ARGUMENTS:
      x, y: numpy ndarrays with (hopefully) matching sizes
            giving coordinates in the old system
      xcen: old-system x coordinate of the new origin
      ycen: old-system y coordinate of the new origin
      phi: angle c.c.w. in degrees from old x to new x axis

    RETURNS: 2-item tuple containing new x and y coordinate arrays

    WRITTEN: Adam S. Bolton, U. of Utah, 2009
    """
    phirad = N.deg2rad(phi)
    xnew = (x - xcen) * N.cos(phirad) + (y - ycen) * N.sin(phirad)
    ynew = (y - ycen) * N.cos(phirad) - (x - xcen) * N.sin(phirad)
    return (xnew,ynew)

def delta_2d(x, y, par):
    """
    NAME: delta_2d

    PURPOSE: Implement 2D Delta function

    USAGE: z = delta_2d(x, y, par)

    ARGUMENTS:
      x, y: vecors or images of coordinates;
            should be matching numpy ndarrays
      par: vector of parameters, defined as follows:
        par[0]: x-center
        par[1]: y-center
        par[2]: radius of a single pixel
        par[3-5]: unused
        
    RETURNS: 2D Delta function evaluated at x-y coords

    WRITTEN: R. S. Bussmann, 2013 November, Cornell University
    """
    square = x.copy()
    square[:] = 0.
    offx = x - par[0]
    offy = y - par[1]
    offset = N.sqrt(offx ** 2 + offy ** 2)
    index = offset < par[2]
    #xspot = x[:, 0] == par[0]
    #yspot = y[0, :] == par[1]
    square[index] = 1.0
    #import pdb; pdb.set_trace()
    #import matplotlib.pyplot as plt
    #plt.imshow(square, origin='lower')
    #plt.show()
    return square

def ellipse_2d(x, y, par):
    """
    NAME: ellipse_2d

    PURPOSE: Implement 2D Elliptical function

    USAGE: z = ellipse_2d(x, y, par)

    ARGUMENTS:
      x, y: vecors or images of coordinates;
            should be matching numpy ndarrays
      par: vector of parameters, defined as follows:
        par[0]: amplitude
        par[1]: intermediate-axis sigma
        par[2]: x-center
        par[3]: y-center
        par[4]: axis ratio
        par[5]: c.c.w. major-axis rotation w.r.t. x-axis
        
    RETURNS: 2D Gaussian evaluated at x-y coords

    NOTE: amplitude = 1 is peak flux, not normalized total flux

    WRITTEN: Adam S. Bolton, U. of Utah, 2009
    """
    (xnew,ynew) = xy_rotate(x, y, -par[2], par[3], par[5])
    r_ell_sq = ((xnew**2)*par[4] + (ynew**2)/par[4]) / N.abs(par[1])**2
    ellipse = r_ell_sq.copy()
    ellipse[:] = 0.
    inside = r_ell_sq < 1
    ellipse[inside] = par[0]
    #import matplotlib.pyplot as plt
    #plt.imshow(r_ell_sq, origin='lower', vmax=10*par[1])
    #plt.colorbar()
    #plt.contour(ellipse)
    #plt.show()
    return ellipse

def gauss_2d(x, y, par):
    """
    NAME: gauss_2d

    PURPOSE: Implement 2D Gaussian function

    USAGE: z = gauss_2d(x, y, par)

    ARGUMENTS:
      x, y: vecors or images of coordinates;
            should be matching numpy ndarrays
      par: vector of parameters, defined as follows:
        par[0]: total flux # in ASB original code, this was amplitude
        par[1]: intermediate-axis sigma
        par[2]: x-center
        par[3]: y-center
        par[4]: axis ratio
        par[5]: c.c.w. major-axis rotation w.r.t. y-axis
        
    RETURNS: 2D Gaussian evaluated at x-y coords

    NOTE: amplitude = 1 is peak flux, not normalized total flux

    WRITTEN: Adam S. Bolton, U. of Utah, 2009
    """
    (xnew,ynew) = xy_rotate(x, y, -par[2], par[3], par[5] + 90)
    r_ell_sq = ((xnew**2)*par[4] + (ynew**2)/par[4]) / N.abs(par[1])**2
    expgauss = N.exp(-0.5*r_ell_sq)
    return par[0] * expgauss

def sie_grad(x, y, par):
    """
    NAME: sie_grad

    PURPOSE: compute the deflection of an SIE potential

    USAGE: (xg, yg) = sie_grad(x, y, par)

    ARGUMENTS:
      x, y: vectors or images of coordinates;
            should be matching numpy ndarrays
      par: vector of parameters with 1 to 5 elements, defined as follows:
        par[0]: lens strength, or 'Einstein radius'
        par[1]: (optional) x-center (default = 0.0)
        par[2]: (optional) y-center (default = 0.0)
        par[3]: (optional) axis ratio (default=1.0)
        par[4]: (optional) major axis Position Angle
                in degrees c.c.w. of y axis. (default = 0.0)

    RETURNS: tuple (xg, yg) of gradients at the positions (x, y)

    NOTES: This routine implements an 'intermediate-axis' convention.
      Analytic forms for the SIE potential can be found in:
        Kassiola & Kovner 1993, ApJ, 417, 450
        Kormann et al. 1994, A&A, 284, 285
        Keeton & Kochanek 1998, ApJ, 495, 157
      The parameter-order convention in this routine differs from that
      of a previous IDL routine of the same name by ASB.

    WRITTEN: Adam S. Bolton, U of Utah, 2009
    """
    # Set parameters:
    b = N.abs(par[0]) # can't be negative!!!
    xzero = 0. if (len(par) < 2) else -par[1]
    yzero = 0. if (len(par) < 3) else par[2]
    q = 1. if (len(par) < 4) else N.abs(par[3])
    phiq = 0. if (len(par) < 5) else par[4]
    eps = 0.001 # for sqrt(1/q - q) < eps, a limit expression is used.
    # Handle q > 1 gracefully:
    if (q > 1.):
        q = 1.0 / q
        phiq = phiq + 90.0
    # Go into shifted coordinats of the potential:
    phirad = N.deg2rad(phiq + 90)
    xsie = (x-xzero) * N.cos(phirad) + (y-yzero) * N.sin(phirad)
    ysie = (y-yzero) * N.cos(phirad) - (x-xzero) * N.sin(phirad)
    # Compute potential gradient in the transformed system:
    r_ell = N.sqrt(q * xsie**2 + ysie**2 / q)
    qfact = N.sqrt(1./q - q)
    # (r_ell == 0) terms prevent divide-by-zero problems
    if (qfact >= eps):
        xtg = (b/qfact) * N.arctan(qfact * xsie / (r_ell + (r_ell == 0)))
        ytg = (b/qfact) * N.arctanh(qfact * ysie / (r_ell + (r_ell == 0)))
        # force r_ell to be greater than 0.1
        thresh = 1e-3
        toolow = (r_ell - b < thresh) & (r_ell - b > 0)
        r_ell[toolow] = r_ell[toolow].mean()#b + thresh
        toolow = (r_ell - b > -thresh) & (r_ell - b < 0)
        r_ell[toolow] = r_ell[toolow].mean()#b - thresh
        M = 1 - b / (r_ell + (r_ell == 0))
        mu = N.abs(1. / M)
        #plt.imshow(mu, origin='lower')
        #plt.colorbar()
        #plt.show()
    else:
        xtg = b * xsie / (r_ell + (r_ell == 0))
        ytg = b * ysie / (r_ell + (r_ell == 0))
        M = 1 - b / (r_ell + (r_ell == 0))
        mu = N.abs(1. / M)
    # Transform back to un-rotated system:
    xg = xtg * N.cos(phirad) - ytg * N.sin(phirad)
    yg = ytg * N.cos(phirad) + xtg * N.sin(phirad)
    #import matplotlib.pyplot as plt
    #plt.imshow(M, origin='lower')
    #plt.colorbar()
    #plt.show()
    # Return value:
    return (xg, yg, mu)

def sbmap(x, y, nlens, nsource, parameters, model_types):

    # define the x, y, and magnification maps
    dx = x.copy()
    dy = y.copy()
    dmu = N.zeros(x.shape)
    nparlens = 5

    # loop over each lens
    for i in range(nlens):

        # Set SIE lens-model parameters and pack them into an array:
        i5 = i * nparlens
        lpar = []
        for ip in range(nparlens):
            lpar.append(parameters[i5 + ip])
        lpar = N.asarray(lpar)

        # Compute the lensing potential gradients and magnification map:
        (xg, yg, mu) = sie_grad(x, y, lpar)

        # apply the gradients and magnifications
        dx -= xg
        dy -= yg
        dmu += mu

    # hack to get the right index from the pzero vector
    interindx = nparlens * nlens

    # loop over each source
    g_lensimage = N.zeros(x.shape)
    g_image = N.zeros(x.shape)
    e_lensimage = N.zeros(x.shape)
    e_image = N.zeros(x.shape)
    amp1 = []
    amp2 = []
    for i in range(nsource):

        i6 = i * 6

        # Set Gaussian source parameters and pack them into an array:
        nparsource = 6
        gpar = []
        for ip in range(nparsource):
            gpar.append(parameters[i6 + interindx + ip])
        gpar = N.asarray(gpar)

        # compute the peak flux of the unlensed Gaussian
        model_type = model_types[i]
        if model_type == 'Delta':
            g_image = delta_2d(x, y, gpar)
        if model_type == 'Gaussian':
            g_image = gauss_2d(x, y, gpar)
        if model_type == 'cylinder':
            g_image = ellipse_2d(x, y, gpar)
        totalflux = g_image.sum()
        if totalflux == 0:
            totalflux = 1.
        normflux = parameters[i6 + interindx] / totalflux
        if model_types[i] != 'Delta':
            gpar[0] *= normflux * 1e-3

        # re-evaluate unlensed image with normalized flux
        if model_type == 'Gaussian':
            g_image = gauss_2d(x, y, gpar)
        if model_type == 'cylinder':
            g_image = ellipse_2d(x, y, gpar)

        if nlens > 0:
            # Evaluate lensed Gaussian image:
            if model_type == 'Delta':
                tmplens = delta_2d(dx, dy, gpar)
            if model_type == 'Gaussian':
                tmplens = gauss_2d(dx, dy, gpar)
            if model_type == 'cylinder':
                tmplens = ellipse_2d(dx, dy, gpar)
            g_lensimage += tmplens
        else:
            # Use the unlensed (but normalized) Gaussian image
            if model_type == 'Delta':
                tmplens = delta_2d(x, y, gpar)
            if model_type == 'Gaussian':
                tmplens = gauss_2d(x, y, gpar)
            if model_type == 'cylinder':
                tmplens = ellipse_2d(x, y, gpar)
            g_lensimage += tmplens

        if nlens > 0:
            # Set elliptical source parameters and pack them into an array:
            epar = gpar.copy()
            epar[0] = 1.0
            epar[1] *= 2.5

            # Evaluate lensed and unlensed elliptical masks:
            lensellipse = ellipse_2d(dx, dy, epar)
            e_lensimage += lensellipse
            ellipse = ellipse_2d(x, y, epar)
            e_image += ellipse

            # Evaluate amplification for each source
            lensmask = lensellipse == 1
            mask = ellipse == 1
            numer = tmplens[lensmask].sum()
            denom = g_image[mask].sum()
            if denom > 0:
                amp_mask = numer / denom
            else:
                amp_mask = 1e2
            numer = tmplens.sum()
            denom = g_image.sum()
            if denom > 0:
                amp_tot = numer / denom
            else:
                amp_tot = 1e2
            if amp_tot > 1e2:
                amp_tot = 1e2
            if amp_mask > 1e2:
                amp_mask = 1e2
            amp1.extend([amp_tot])
            amp2.extend([amp_mask])

    return g_image, g_lensimage, e_image, e_lensimage, amp1, amp2
