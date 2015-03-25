# encoding: utf-8
# cython: profile=True
# filename: sample_vis.pyx
"""
2013 March 29
Shane Bussmann
Ported from an IDL routine of the same name from Katherine Rosenfeld
 - Sampling the model visibilities using the UV of the observation
 - Trying to implement Schwab's "Optimal Gridding" system that MIRIAD uses
    F.R. Schwab "Optimal gridding of visibility data in radio interferometry",
      Indirect Imaging, J.A. Roberts, ed., Cambridge Univ. Press, 1984.
 - See MIRIAD's ModMap in model.for
"""

from __future__ import print_function, division

# load constants and helper functions
#@grid.pro
#@ModGrid.pro
#@constants.pro
#import pyfits
#import math
#import ModGrid
import numpy as np
import grid
#cimport np

# define some constants
c = 3e10
kB = 1.3808e-16
pc = 3.08572e18
Jy = 1e23
pi = np.pi
Lsun = 3.826e33
Msun = 2e33
G = 6.67e-8
rad = 206265.
kms = 1e5
GHz = 1e9


#DTYPE = np.double
#ctypedef np.double_t DTYPE_t

def ModGrid1(vv, uu, v0, u0, Grd, gcf, ngcf, nyd, width):

    #assert Grd.dtype == DTYPE and gcf.dtype == DTYPE# and uuu.dtype == DTYPE and vvv.dtype == DTYPE

    #cdef:
        #int ju, jv, p, q, jun2, jun1, jvn2, jvn1, jup1, jup2, jup3, i
        #int p, q
        #double u, v, wu1, wu2, wu3, wu4, wu5, wu6, wv1, wv2, wv3, wv4, wv5, \
        #        wv6, w, Intpi
        #np.ndarray Intp

    uuu = uu.copy()
    vvv = vv.copy()
    #conjugate = (uuu >= 0)
    #uuu[conjugate] = u0 + 1 * uuu[conjugate]
    #vvv[conjugate] = v0 + 1 * vvv[conjugate]
    conjugate = (uuu < 0)
    uuu[conjugate] = -1 * uuu[conjugate]
    vvv[conjugate] = -1 * vvv[conjugate]
    #vvv[conjugate] = nyd - np.abs(vvv[conjugate])
    negv = (vvv < 0)
    bugstep = np.abs(vvv[negv])
    vvv[negv] = nyd - bugstep

    # np bug(?): vectors with values close to zero are handled poorly
    #checkhigh = vvv == nyd
    #vvv[checkhigh] = nyd - 1
    #print(uuu)
    #print(vvv)

    nvis = vvv.size
    wu = np.zeros([width, nvis])
    wv = np.zeros([width, nvis])
    step = (ngcf - 1) // width
    rv = vvv - np.floor(vvv)
    ru = uuu - np.floor(uuu)
    uuu = np.floor(uuu).astype(int)
    vvv = np.floor(vvv).astype(int)
    p = ngcf // 2 - np.around(step * rv)
    q = ngcf // 2 - np.around(step * ru)
    #print(p)
    #print(q)
    #toohigh = p >= ngcf // 2
    #p[toohigh] = p[toohigh] - ngcf
    #toohigh = q >= ngcf / 2
    #q[toohigh] = q[toohigh] - ngcf
    for i in range(width):
        wuindx = q.astype(int) + step * (i - width // 2 + 1)
        wu[i, :] = gcf[wuindx]
        wvindx = p.astype(int) + step * (i - width // 2 + 1)
        wv[i, :] = gcf[wvindx]
        #print(uuu.min(), uuu.max(), vvv.min(), vvv.max())
        #print(wuindx - ngcf // 2, gcf[wuindx])

    #print(wu.sum(axis=0).shape, wv.sum(axis=0).shape)
    w = wu.sum(axis = 0) * wv.sum(axis = 0)
    #print(w.shape)
    #print(uuu.max(), vvv.max())

    #uuun2 = uuu - 2
    #uuun1 = uuu - 1
    #uuup1 = uuu + 1
    #uuup2 = uuu + 2
    #uuup3 = uuu + 3
    #vvvn2 = vvv - 2
    #vvvn1 = vvv - 1
    #vvvp1 = vvv + 1
    #vvvp2 = vvv + 2
    #vvvp3 = vvv + 3

    # if nyd is 128, then if vvv is 127, it becomes -1.  vvv=125 --> -3
    toohigh = vvv >= nyd - 3
    vvv[toohigh] = vvv[toohigh] - nyd
    #toohigh = vvvp1 >= nyd
    #vvvp1[toohigh] = vvvp1[toohigh] - nyd
    #toohigh = vvvp2 >= nyd
    #vvvp2[toohigh] = vvvp2[toohigh] - nyd
    #toohigh = vvvp3 >= nyd
    #vvvp3[toohigh] = vvvp3[toohigh] - nyd
    #toolow = uuun2 < 0
    #uuun2[toolow] = nyd + (uuun2[toolow])
    #toolow = uuun1 < 0
    #uuun1[toolow] = nyd + (uuun1[toolow])
    #toolow = vvvn2 < 0
    #vvvn2[toolow] = nyd + (vvvn2[toolow])
    #toolow = vvvn1 < 0
    #vvvn1[toolow] = nyd + (vvvn1[toolow])
    #    if ju - 2 < 0:
    #        jun2 = nyd + (ju - 2)
    #    else:
    #        jun2 = ju - 2
    #    if ju - 1 < 0:
    #        jun1 = nyd + (ju - 1)
    #    else:
    #        jun1 = ju - 1
    #    if jv - 2 < 0:
    #        jvn2 = nyd + (jv - 2)
    #    else:
    #        jvn2 = jv - 2
    #    if jv - 1 < 0:
    #        jvn1 = nyd + (jv - 1)
    #    else:
    #        jvn1 = jv - 1

        #if jv + 1 >= nyd:
        #    jvp1 = jv + 1 - nyd
        #else:
        #    jvp1 = jv + 1 
        #if jv+2 >= nyd:
        #    jvp2 = jv + 2 - nyd
        #else:
        #    jvp2 = jv + 2 
        #if jv + 3 >= nyd:
        #    jvp3 = jv + 3 - nyd
        #else:
        #    jvp3 = jv + 3

    Intp = np.zeros(nvis, dtype='complex')
    for j in range(width):
        ivvv = vvv - (width - 1) // 2 + j
        for i in range(width):
            iuuu = uuu - (width - 1) // 2 + i
            Intp[:] += 1 / w * wv[j, :] * wu[i, :] * Grd[iuuu, ivvv]

    #  Conjugate the data if necessary.
    Intp[conjugate] = np.conjugate(Intp[conjugate])

    return Intp


def uvmodel(model, modelheader, u, v, pcd):

    #model = ''
    #modelheader = ''
    #u = ''
    #v = ''
    #visheader = ''
    # - set up information
    # - files
    #imfile = 'g_lensimage.fits'		# - model image file
    #imvfile = 'g_lensimage_ext_day1_lsb.uvfits'		# - visibilities produce by MIRIAD
    #visfile = 'G09v1.97_ext_day1_lsb.uvfits'	# - data visibility file (CONTINUUM)

     # - grid settings
    igrid = 1	 # - (1) weighted (2) nearest neighbor
    alpha = 1.	 # - hard coded
    width = 6	 # - hard coded
    maxgcf = 2048

    # - griding info: read in data visibilities
    #miriad = pyfits.getdata(visfile)
    #visheader = pyfits.getheader(visfile)
    #freq = visheader['crval4']

    # - read in model file
    #model = pyfits.getdata(imfile)	# transpose to mimic how array is read in
    #modelheader = pyfits.getheader(imfile)

    # - FFT parameters (follow ModMap variables)
    nx = modelheader['NAXIS1']	# - # of x pixels
    ny = modelheader['NAXIS2']	# - # of y pixels
    ln2 = np.log(2)
    nxd = 2 ** np.ceil(np.log(2. * nx) / ln2)	# - padding w/ zeros
    nyd = 2 ** np.ceil(np.log(2. * ny) / ln2)	# - padding w/ zeros
    nxd = np.long(nxd)
    nyd = np.long(nyd)
    dx = modelheader['CDELT1'] * pi / 180  # - pixel size in radians
    dy = modelheader['CDELT2'] * pi / 180  # - pixel size in radians 
    du = 1. / (dx * nxd)		# - size of u cell
    dv = 1. / (dy * nyd)		# - size of v cell
    nu = nxd			# - # of u cells
    nv = nyd			# - # of v cells
    u0 = nu / 2 + 1
    v0 = nv / 2 + 1
    """
    Is this necessary?!?
    umax = 0.5 * (nxd - 1 - width)
    vmax = 0.5 * (nyd - 1 - width)
    """

    # - now follow ModFFT in model.for
    #dra  = -1e0 * modelheader['CDELT1'] * pi / 180  # - pixel size in radians 
    #ddec = modelheader['CDELT2'] * pi / 180  # - pixel size in radians 
    iref = nx // 2 + 1
    jref = ny // 2 + 1
    modelcrpix1 = modelheader['CRPIX1']
    modelcrpix2 = modelheader['CRPIX2']
    pcm = [modelheader['CRVAL1'], modelheader['CRVAL2']]
    raref1 = dx * (modelcrpix1 - iref)
    decref1 = dy * (modelcrpix2 - jref)
    raref2 = (pcd[0] - pcm[0]) * \
        np.cos(pcd[1] * pi / 180.) * pi / 180.
    decref2 = (pcd[1] - pcm[1]) * pi / 180.

    # - CALCULATE UVW CONVERSION MATRIX FOR GEOMETRY DIFFERENCES
    # phase center for the model
    ucoeff, vcoeff, wcoeff = grid.coGeom(pcm, pcd)

    # - shift array and pad with zeros
    shifted = np.zeros([nyd, nxd])
    noff = nxd - nx
    areal = shifted.copy()
    aimag = shifted.copy()
    #x1 = nxd // 2 - nx // 2
    #x2 = nxd // 2 + nx // 2
    #y1 = nyd // 2 - ny // 2
    #y2 = nyd // 2 + ny // 2
    #areal[x1:x2, y1:y2] = model
    areal[noff:, noff:] = model
    #image = model + 0j * model
    image = areal + 1j * aimag
    image = np.roll(image, ny // 2, axis=0)
    image = np.roll(image, nx // 2, axis=1)

    #mu_grid, dump = np.meshgrid(mu, np.zeros(nxd) + 1)
    #mu_grid_shifty = np.roll(mu_grid, -1 * u0int, axis=0)
    #mu_grid_shiftyx = np.roll(mu_grid, -1 * u0int, axis=1)
    #foo = du * mu_grid_shiftyx
    #tmu = np.transpose(mu)
    #tmu_grid, dump = np.meshgrid(tmu, np.zeros(nxd) + 1)
    #tmu_grid_shifty = np.roll(tmu_grid, -1 * u0int, axis=0)
    #tmu_grid_shiftyx = np.roll(tmu_grid, -1 * u0int, axis=1)
    #mv = du * tmu_grid_shiftyx
    #mu = foo[::-1,:]

    # - read in the data visibilities
    #u = freq * miriad['UU']
    #v = freq * miriad['VV']
    ud = ucoeff[0] * u + ucoeff[1] * v
    vd = vcoeff[0] * u + vcoeff[1] * v

    #print(u.max(), ud.max(), du)

    nvis = len(u)
    #print('# of visibilites:', nvis)

    # - calculate gridding correction function
    if igrid == 1:
        #print('ModCorr (grid.for)')
        #import time
        #start = time.time()
        ycorr, xcorr = grid.ModCorr(nyd, nxd)
        mcorr = np.outer(ycorr, xcorr)
        image = image / mcorr
        #time_modgrid = time.time()-start
        #print('time to run ModCorr: ', time_modgrid, 'seconds')

    # - take FFT
    #print('ModPlane (model.for)')
    #cdef np.ndarray mvis_real
    #cdef np.ndarray mvis_imag
    mvis = np.fft.fft2(image)#, s=(nyd, nxd))
    mvis = np.conjugate(mvis)

    # - Follow ModGrid in model.for
    #print('ModGrid (model.for)')
    #cdef np.ndarray gcf
    ngcf = width * ((maxgcf - 1) // width) + 1
    gcf = grid.gcffun(ngcf, width, alpha)

    #cdef np.ndarray uu
    #cdef np.ndarray vv

    # - interpolate grid
    if igrid == 1:
        #print('interpolating grid')
        # mvis = transpose(mvis)
        uu = ud / du  # calculate fractional index in the grid
        vv = vd / dv  # calculate fractional index in the grid
        #mvis_opt = np.zeros(nvis) + 1.j * np.zeros(nvis)
        # - weighted interpolation
        #import pyximport
        #pyximport.install()
        #import ModGrid
        #for i in np.arange(nvis):
        #    if (np.abs(uu[i]) > umax) or (np.abs(vv[i]) > vmax):
        #        mvis_opt[i] = 0.
        #    else:
        #        start = time.time()
        #mvis_real = np.real(mvis)
        #mvis_imag = np.imag(mvis)
        #uuint = np.around(uu).astype(int)
        #vvint = np.around(vv).astype(int)
        #ru = uu - uuint
        #rv = vv - vvint
        mvis_opt = ModGrid1(uu, vv, u0, v0, mvis, gcf, ngcf, nyd, width)
        #print(uu)
        #print(vv)
        """
        Is this necessary?!?
        overmax = np.abs(vv) > vmax
        mvis_opt[overmax] = 0.
        overmax = np.abs(uu) > umax
        mvis_opt[overmax] = 0.
        """
        #mvis_opt1 = np.zeros(nvis) + 1.j * np.zeros(nvis)
        #p = np.zeros(nvis)
        #q = np.zeros(nvis)
        #for i in np.arange(nvis):
        #    if (np.abs(uu[i]) > umax) or (np.abs(vv[i]) > vmax):
        #        mvis_opt1[i] = 0.
        #    else:
        #        start = time.time()
        #        mvis_real = np.real(mvis)
        #        mvis_imag = np.imag(mvis)
        #        mvis_opt1[i], p[i], q[i] = ModGrid.ModGrid(uu[i], vv[i], mvis, 1, 1, 1, \
        #                gcf, ngcf, nyd)
        #        time_modgrid = time_modgrid + time.time()-start
        #        start = time.time()
                #mvis_opt[i] = grid.ModShift(ud[i], vd[i], raref1, decref1, raref2, \
                #        decref2, 1, 1, mvis_opt[i])

        #print(uu)
        #diffmvis = np.abs((mvis_opt - mvis_opt1) / mvis_opt)
        #print(diffmvis)
        #print(diffmvis.min(), diffmvis.max())
        #mvis_opt_imag = ModGrid(uuint, vvint, mvis_imag, gcf, ngcf, nyd, width, step)
        #mvis_opt = mvis_opt_real + 1.j * mvis_opt_imag
        #        time_modgrid = time_modgrid + time.time()-start
        #        start = time.time()
        mvis_opt = grid.ModShift(ud, vd, raref1, decref1, raref2, \
                        decref2, 1, 1, mvis_opt)
        #for i in np.arange(nvis):
        #    if (np.abs(uu[i]) > umax) or (np.abs(vv[i]) > vmax):
        #        mvis_opt[i] = 0.
        #    else:
        #        mvis_opt[i] = grid.ModShift(ud[i], vd[i], raref1, decref1, raref2, \
        #                decref2, 1, 1, mvis_opt[i])
                #time_modshift = time_modshift + time.time()-start
        #print(time_modgrid, time_modshift)
    if igrid == 2:
        u0 = nxd // 2 + 1		# - central u cell id
        #v0 = nyd // 2 + 1		# - central v cell id
        mu = np.arange(nxd)
        u0int = np.int(u0)
        mu[u0int:] = u0int - nxd + np.arange(u0int - 2)
        mu1d = du * mu
        import time
        # - Nearest neighbor search
        uu = ud / du #  calculate index in the grid
        vv = vd / dv #  calculate index in the grid
        mvis_opt = np.zeros(nvis) + 1.j * np.zeros(nvis)
        mv1d = -1 * mu1d
        imu1d = np.argsort(mu1d)
        smu1d = mu1d[imu1d]
        imv1d = np.argsort(mv1d)
        smv1d = mv1d[imv1d]
        time_modgrid = 0
        for i in np.arange(nvis):
            start = time.time()
            iu = (np.abs(smu1d - u[i])).argmin()
            iv = (np.abs(smv1d - v[i])).argmin()
            mvis_opt[i] = mvis[imv1d[iv], imu1d[iu]]
            mvis_opt[i] = grid.ModShift(ud[i], vd[i], raref1, decref1, raref2, \
                    decref2, 1, 1, mvis_opt[i])
            time_modgrid = time_modgrid + time.time()-start
        #print(time_modgrid)

    #mvis_opt = np.conjugate(mvis_opt)
    # - read in visibilities MIRIAD produces
    #imvis = pyfits.getdata(imvfile)
    #mirvisheader = pyfits.getheader(imvfile)
    #ruv = np.sqrt(u ** 2 + v ** 2)

    # read c version
    #readcol,'../myvis.dat',cu,cv,cr,ci

    # - plot the results
    #rvis = imvis['DATA'][:, 0, 0, 0, 0, 0]
    #ivis = imvis['DATA'][:, 0, 0, 0, 0, 1]
    #wgt = imvis['DATA'][:, 0, 0, 0, 0, 2]
    #ampvis = np.sqrt(rvis ** 2 + ivis ** 2)

    # Real part
    #mvis_amp = np.sqrt(np.real(mvis_opt) ** 2 + np.imag(mvis_opt) ** 2)
    #mvis_amp = np.abs(mvis_opt)
    #plt.clf()
    #plt.plot(ruv, mvis_amp-ampvis, ',', color='red', zorder=2)
    #plt.axis([0, 300, 0, 1])
    #plt.plot(ruv, np.real(mvis_opt), 'o', ms=2, mec='red', mfc='red')  # Python
    #plt.plot(ruv, rvis, '.', ms=1, color='black')  # MIRIAD
    #plt.show()
    #saveloc = 'real_python_miriad'
    #savefig(saveloc)

     # Imaginary part
    #import matplotlib.pyplot as plt
    #from pylab import savefig
    #plt.clf()
    #plt.plot(ruv, np.abs(mvis_opt), 'o', ms=2, mec='red', mfc='red')  # Python
    #plt.loglog()
    ##plt.plot(ruv, ivis, 'o', mfc='black', ms=1)  # MIRIAD
    #saveloc = 'imaginary_python_miriad'
    #savefig(saveloc)

    # calculate chi-squared 
    #chi2v = wgt * (np.real(mvis_opt) - rvis) ** 2 + \
    #    wgt * (np.imag(mvis_opt) - ivis) ** 2
    #print('chi squared value is :', chi2v.sum())
    # 4.7212017638007328e-06
    return mvis_opt
