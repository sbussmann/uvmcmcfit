# - return gridding correction for Mr. Schwab's optimal gridding system
# - Essentially a transription of grid.for

from __future__ import print_function

import numpy
import math

def getspherewave():

    # Define the spheroidal weighting function grid (See MIRIAD's grid.for)
    # Calculates the spheroidal wave function.
    # These rational approximations are taken from Schwab "Optimal Gridding",
    #  Indirect Imaging, ed J.A. Robert, 1984. 

    p = numpy.zeros([7, 5, 5, 2])
    q = numpy.zeros([3, 5, 5, 2])

    # M=4, ALPHA=0,2,0.4, ETA<ETALIM
    p[:,:,0,0] =   numpy.array( \
            [1.584774E-2,-1.269612E-1, 2.333851E-1, \
            -1.636744E-1, 5.014648E-2, 0.0, 0.0, \
             3.101855E-2,-1.641253E-1, 2.385500E-1, \
            -1.417069E-1, 3.773226E-2, 0.0, 0.0, \
             5.007900E-2,-1.971357E-1, 2.363775E-1, \
            -1.215569E-1, 2.853104E-2, 0.0, 0.0, \
             7.201260E-2,-2.251580E-1, 2.293715E-1, \
            -1.038359E-1, 2.174211E-2, 0.0, 0.0, \
             9.585932E-2,-2.481381E-1, 2.194469E-1, \
            -8.862132E-2, 1.672243E-2, 0.0, 0.0]).reshape(7, 5, order='F')
    # - M=5, ALPHA=0,2,0.5, ETA<ETALIM
    p[:,:,1,0]  =   numpy.array( \
            [3.722238E-3,-4.991683E-2, 1.658905E-1,-2.387240E-1, \
             1.877469E-1,-8.159855E-2, 3.051959E-2, \
             8.182649E-3,-7.325459E-2, 1.945697E-1,-2.396387E-1, \
             1.667832E-1,-6.620786E-2, 2.224041E-2, \
             1.466325E-2,-9.858686E-2, 2.180684E-1,-2.347118E-1, \
             1.464354E-1,-5.350728E-2, 1.624782E-2, \
             2.314317E-2,-1.246383E-1, 2.362036E-1,-2.257366E-1, \
             1.275895E-1,-4.317874E-2, 1.193168E-2, \
             3.346886E-2,-1.503778E-1, 2.492826E-1,-2.142055E-1, \
             1.106482E-1,-3.486024E-2, 8.821107E-3]).reshape(7, 5, order='F')
    # M=6, ALPHA=0,2,0.5, ETA<ETALIM
    p[:,:,2,0]  =   numpy.array( \
            [5.613913E-2,-3.019847E-1, 6.256387E-1, \
            -6.324887E-1, 3.303194E-1, 0.0, 0.0, \
             6.843713E-2,-3.342119E-1, 6.302307E-1, \
            -5.829747E-1, 2.765700E-1, 0.0, 0.0, \
             8.203343E-2,-3.644705E-1, 6.278660E-1, \
                -5.335581E-1, 2.312756E-1, 0.0, 0.0, \
                 9.675562E-2,-3.922489E-1, 6.197133E-1, \
                -4.857470E-1, 1.934013E-1, 0.0, 0.0, \
                 1.124069E-1,-4.172349E-1, 6.069622E-1, \
                -4.405326E-1, 1.618978E-1, 0.0, 0.0]).reshape(7, 5, order='F')
    # M=7, ALPHA=0,2,0.5, ETA<ETALIM
    p[:,:,3,0]  =   numpy.array( \
            [2.460495e-2,-1.640964e-1, 4.340110e-1, \
            -5.705516e-1, 4.418614e-1, 0.0, 0.0, \
             3.070261e-2,-1.879546e-1, 4.565902e-1, \
            -5.544891e-1, 3.892790e-1, 0.0, 0.0, \
             3.770526e-2,-2.121608e-1, 4.746423E-1, \
            -5.338058e-1, 3.417026e-1, 0.0, 0.0, \
             4.559398e-2,-2.362670e-1, 4.881998e-1, \
            -5.098448e-1, 2.991635e-1, 0.0, 0.0, \
             5.432500e-2,-2.598752e-1, 4.974791e-1, \
            -4.837861e-1, 2.614838e-1, 0.0, 0.0]).reshape(7, 5, order='F')
    # M=8, ALPHA=0,2,0.5, ETA<ETALIM				
    p[:,:,4,0]  =   numpy.array( \
            [1.378030e-2,-1.097846e-1, 3.625283e-1, \
            -6.522477e-1, 6.684458e-1,-4.703556e-1,0.0, \
             1.721632e-2,-1.274981e-1, 3.917226e-1, \
            -6.562264e-1, 6.305859e-1,-4.067119e-1,0.0, \
             2.121871e-2,-1.461891e-1, 4.185427e-1, \
            -6.543539e-1, 5.904660e-1,-3.507098e-1,0.0, \
             2.580565e-2,-1.656048e-1, 4.426283e-1, \
            -6.473472e-1, 5.494752e-1,-3.018936e-1,0.0, \
             3.098251e-2,-1.854823e-1, 4.637398e-1, \
            -6.359482e-1, 5.086794e-1,-2.595588e-1,0.0]).reshape(7, 5, order='F')
    # M=6, ALPHA=0,2,0.5, ETA>ETALIM
    p[:,:,2,1]  =   numpy.array( \
            [8.531865E-4,-1.616105E-2, 6.888533E-2, \
                -1.109391E-1, 7.747182E-2, 0.0, 0.0, \
                 2.060760E-3,-2.558954E-2, 8.595213E-2, \
                -1.170228E-1, 7.094106E-2, 0.0, 0.0, \
                 4.028559E-3,-3.697768E-2, 1.021332E-1, \
                -1.201436E-1, 6.412774E-2, 0.0, 0.0, \
                 6.887946E-3,-4.994202E-2, 1.168451E-1, \
                -1.207733E-1, 5.744210E-2, 0.0, 0.0, \
                 1.071895E-2,-6.404749E-2, 1.297386E-1, \
                -1.194208E-1, 5.112822E-2, 0.0, 0.0]).reshape(7, 5, order='F')
    # M=7, ALPHA=0,2,0.5, ETA>ETALIM
    p[:,:,3,1]  =   numpy.array( \
            [1.924318e-4,-5.044864e-3, 2.979803e-2, \
            -6.660688e-2, 6.792268e-2, 0.0, 0.0, \
             5.030909e-4,-8.639332e-3, 4.018472e-2, \
            -7.595456e-2, 6.696215e-2, 0.0, 0.0, \
             1.059406e-3,-1.343605e-2, 5.135360e-2, \
            -8.386588e-2, 6.484517e-2, 0.0, 0.0, \
             1.941904e-3,-1.943727e-2, 6.288221e-2, \
            -9.021607e-2, 6.193000e-2, 0.0, 0.0, \
             3.224785e-3,-2.657664e-2, 7.438627e-2, \
            -9.500554e-2, 5.850884e-2, 0.0, 0.0]).reshape(7, 5, order='F')
    # M=8, ALPHA=0,2,0.5, ETA>ETALIM
    p[:,:,4,1]   =   numpy.array( \
            [4.290460e-5,-1.508077e-3, 1.233763e-2, \
            -4.091270e-2, 6.547454e-2,-5.664203e-2,0.0, \
             1.201008e-4,-2.778372e-3, 1.797999e-2, \
            -5.055048e-2, 7.125083e-2,-5.469912e-2,0.0, \
             2.698511e-4,-4.628815e-3, 2.470890e-2, \
            -6.017759e-2, 7.566434e-2,-5.202678e-2,0.0, \
             5.259595e-4,-7.144198e-3, 3.238633e-2, \
            -6.946769e-2, 7.873067e-2,-4.889490e-2,0.0, \
             9.255826e-4,-1.038126e-2, 4.083176e-2, \
            -7.815954e-2, 8.054087e-2,-4.552077e-2,0.0]).reshape(7, 5, order='F')
    # M=4, ALPHA=0,2,0.5, ETA<ETALIM
    q[:,:,0,0]   = numpy.array( \
            [1., 4.845581E-1, 7.457381E-2, \
            1., 4.514531E-1, 6.458640E-2, \
            1., 4.228767E-1, 5.655715E-2, \
            1., 3.978515E-1, 4.997164E-2, \
            1., 3.756999E-1, 4.448800E-2]).reshape(3, 5, order='F')
    # M=5, ALPHA=0,2,0.5, ETA<ETALIM
    q[:,:,1,0]   = numpy.array( \
            [1., 2.418820E-1, 0.0, \
            1., 2.291233E-1, 0.0, \
            1., 2.177793E-1, 0.0, \
            1., 2.075784E-1, 0.0, \
            1., 1.983358E-1, 0.0]).reshape(3, 5, order='F')
    # M=6, ALPHA=0,2,0.5, ETA<ETALIM
    q[:,:,2,0]   = numpy.array( \
            [1., 9.077644E-1, 2.535284E-1, \
            1., 8.626056E-1, 2.291400E-1, \
            1., 8.212018E-1, 2.078043E-1, \
            1., 7.831755E-1, 1.890848E-1, \
            1., 7.481828E-1, 1.726085E-1]).reshape(3, 5, order='F')
    # M=7, ALPHA=0,2,0.5, ETA<ETALIM
    q[:,:,3,0]   = numpy.array( \
            [1., 1.124957e00, 3.784976e-1, \
            1., 1.075420e00, 3.466086e-1, \
            1., 1.029374e00, 3.181219e-1, \
            1., 9.865496e-1, 2.926441e-1, \
            1., 9.466891e-1, 2.698218e-1]).reshape(3, 5, order='F')
    # M=7(8?), ALPHA=0,2,0.5, ETA<ETALIM
    q[:,:,4,0]   = numpy.array( \
            [1., 1.076975e00, 3.394154e-1, \
            1., 1.036132e00, 3.145673e-1, \
            1., 9.978025e-1, 2.920529e-1, \
            1., 9.617584e-1, 2.715949e-1, \
            1., 9.278774e-1, 2.530051e-1]).reshape(3, 5, order='F')
    # M=6, ALPHA=0,2,0.5, ETA>ETALIM
    q[:,:,2,1]   = numpy.array( \
            [1., 1.101270   , 3.858544E-1, \
            1., 1.025431   , 3.337648E-1, \
            1., 9.599102E-1, 2.918724E-1, \
            1., 9.025276E-1, 2.575337E-1, \
            1., 8.517470E-1, 2.289667E-1]).reshape(3, 5, order='F')
    # M=7, ALPHA=0,2,0.5, ETA>ETALIM
    q[:,:,3,1]   = numpy.array( \
            [1., 1.450730e00, 6.578684e-1, \
            1., 1.353872e00, 5.724332e-1, \
            1., 1.269924e00, 5.032139e-1, \
            1., 1.196177e00, 4.460948e-1, \
            1., 1.130719e00, 3.982785e-1]).reshape(3, 5, order='F')
    # M=8, ALPHA=0,2,0.5, ETA>ETALIM
    q[:,:,4,1]   = numpy.array( \
            [1., 1.379457e00, 5.786953e-1, \
            1., 1.300303e00, 5.135748e-1, \
            1., 1.230436e00, 4.593779e-1, \
            1., 1.168075e00, 4.135871e-1, \
            1., 1.111893e00, 3.744076e-1]).reshape(3, 5, order='F')

    return p, q

def spheroid(eta, m, alpha, p, q):
    # Currently alpha can only be eq to 1
    if (alpha != 1):
        print('grid.py: ALPHA MUST BE 1')

    etalim = [1., 1., 0.75, 0.775, 0.775]
    nnum   = [5, 7, 5, 5, 6]
    ndenom = [3, 2, 3, 3, 3]

    # checks and balances
    twoalp = numpy.int(numpy.round(2. * alpha))
    if (numpy.abs(eta) > 1):
        print("Abs(ETA) exceeds 1: {:f}".format(eta))
    if (twoalp < 0) or (twoalp > 4):
        print("Illegal value of ALPHA")
    if (m < 4) or (m > 8):
        print("Illegal value of M: {:f}".format(m))

    # - Go to appropriate approximation
    if (numpy.abs(eta) > etalim[m - 4]):
        ip = 2 - 1
        x = eta * eta - 1
    else:
        ip = 1 - 1
        x  = eta * eta - etalim[m - 4] * etalim[m - 4]

    # - Get numerator via Horners rule:
    np = nnum[m - 4] - 1
    num = p[np, twoalp, m - 4, ip]
    for i in range(np - 1, -1, -1):
        num = num * x + p[i, twoalp, m - 4, ip]

    # - Get denominator via Horners rule"
    nq = ndenom[m - 4] - 1
    denom = q[nq, twoalp, m - 4, ip] 
    for i in range(nq - 1, -1, -1):
        denom = denom * x + q[i, twoalp, m - 4, ip]

    return num/denom

def gcffun(n, width, alpha):
    ppp, qqq = getspherewave()
    phi = numpy.zeros(n)
    j = numpy.int(numpy.round(2. * alpha))
    p = 0.5 * j
    if (j == 0):
        for i in numpy.arange(n):
            x = (2. * i - (n - 1.)) / (n - 1.)
            phi[i] = spheroid(x, width, p, ppp, qqq)
    else:
        for i in numpy.arange(n):
            x = (2. * i - (n - 1.)) / (n - 1.)
            phi[i] = numpy.sqrt(1. - x * x) ** j * spheroid(x, width, p, ppp, qqq)
    return phi

def corrfun(n, width, alpha):
# See MIRIAD's grid.for
    ppp, qqq = getspherewave()
    phi = numpy.zeros(n)
    dx = 2. / n
    i0 = numpy.int(n) // 2 + 1
    for i in range(n):
        x = (i + 1 - i0) * dx
        phi[i] = spheroid(x, width, alpha, ppp, qqq)
    return phi

def ModCorr(nxd, nyd):
    # See MIRIAD's model.for
    #	- which include half-image shift in a (-1)**j-1 factor
    width = 6

    xcorr = numpy.zeros(nxd)
    #xcorr1 = numpy.zeros(nxd)
    ycorr = numpy.zeros(nyd)
    #ycorr1 = numpy.zeros(nyd)

    data = corrfun(nxd, width, 1.)
    offset = numpy.int(nxd) // 2
    indx = numpy.arange(nxd // 2)
    ix = indx.astype(int)
    xcorr[ix] = data[ix + offset]
    indx = numpy.arange(nxd // 2) + nxd // 2
    ix = indx.astype(int)
    xcorr[ix] = data[ix - offset]
    #for i in numpy.arange(nxd // 2):
    #    xcorr1[i] = data[i + offset]
    #for i in numpy.arange(nxd // 2) + nxd // 2:
    #    xcorr1[i] = data[i - offset]

    data = corrfun(nyd, width, 1.)
    offset = numpy.int(nyd) // 2
    indx = numpy.arange(0, nyd // 2, 2)
    ycorr[indx] = data[indx + offset]
    indx = numpy.arange(0, nyd // 2, 2)
    ycorr[indx + 1] = data[indx + offset + 1]
    indx = numpy.arange(nyd // 2, nyd, 2)
    ycorr[indx] = data[indx - offset]
    indx = numpy.arange(nyd // 2, nyd, 2)
    ycorr[indx + 1] = data[indx - offset + 1]
    #for i in numpy.arange(0, nyd // 2, 2):
    #    ycorr1[i]   =  data[i + offset]
    #    ycorr1[i + 1] =  data[i + 1 + offset]
    #for i in numpy.arange(nyd // 2, nyd, 2):
    #    ycorr1[i]   =  data[i - offset]
    #    ycorr1[i + 1] =  data[i + 1 - offset]
    #print(ycorr - ycorr1)
    return ycorr, xcorr

def ModShift(uu, vv, xref1, yref1, xref2, yref2, freq1, freq, intp):
    uu = numpy.array(uu)
    vv = numpy.array(vv)
    t1 = -2. * math.pi * (uu * xref1 + vv * yref1)
    t2 = -2. * math.pi * (uu * xref2 + vv * yref2) / freq1
    theta = t1 + t2 * freq
    #W = numpy.complex(numpy.cos(theta), numpy.sin(theta))
    W = numpy.cos(theta) + 1.j * numpy.sin(theta)
    return W * intp

def coGeom(phase_center_model, phase_center_data):
    # See coGeom in subs/co.for
    #@constants.pro	
    ucoeff = numpy.zeros(3)
    vcoeff = numpy.zeros(3)
    wcoeff = numpy.zeros(3) 

    # calculate model coordinates
    lng0 = phase_center_model[0] * math.pi / 180
    lat0 = phase_center_model[1] * math.pi / 180
    # calculate data coordinates
    lng = phase_center_data[0] * math.pi / 180
    lat = phase_center_data[1] * math.pi / 180
    # compute matrix elements
    clat0 = numpy.cos(lat0)
    slat0 = numpy.sin(lat0)
    clng  = numpy.cos(lng - lng0)
    slng  = numpy.sin(lng - lng0)
    clat  = numpy.cos(lat)
    slat  = numpy.sin(lat)

    # CODE MUST BE TYPE SIN
    fac = 1e0 / (slat * slat0 + clat * clat0 * clng)
    ucoeff[0] =  fac * (clat * clat0 + slat * slat0 * clng)
    ucoeff[1] = -fac * slat0 * slng
    ucoeff[2] =  0.
    vcoeff[0] =  fac * slat * slng
    vcoeff[1] =  fac * clng
    vcoeff[2] =  0.
    wcoeff[0] =  0.
    wcoeff[1] =  0.
    wcoeff[2] =  0.

    return ucoeff, vcoeff, wcoeff
