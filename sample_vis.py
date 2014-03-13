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
import numpy
import grid
#cimport numpy
#import pdb

# define some constants
c = 3e10
kB = 1.3808e-16
pc = 3.08572e18
Jy = 1e23
pi = numpy.pi
Lsun = 3.826e33
Msun = 2e33
G = 6.67e-8
rad = 206265.
kms = 1e5
GHz = 1e9


#DTYPE = numpy.double
#ctypedef numpy.double_t DTYPE_t

def ModGrid1(vv, uu, v0, u0, Grd, gcf, ngcf, nyd, width):

    #assert Grd.dtype == DTYPE and gcf.dtype == DTYPE# and uuu.dtype == DTYPE and vvv.dtype == DTYPE

    #cdef:
        #int ju, jv, p, q, jun2, jun1, jvn2, jvn1, jup1, jup2, jup3, i
        #int p, q
        #double u, v, wu1, wu2, wu3, wu4, wu5, wu6, wv1, wv2, wv3, wv4, wv5, \
        #        wv6, w, Intpi
        #numpy.ndarray Intp

    uuu = uu.copy()
    vvv = vv.copy()
    #conjugate = (uuu >= 0)
    #uuu[conjugate] = u0 + 1 * uuu[conjugate]
    #vvv[conjugate] = v0 + 1 * vvv[conjugate]
    conjugate = (uuu < 0)
    uuu[conjugate] = -1 * uuu[conjugate]
    vvv[conjugate] = -1 * vvv[conjugate]
    #vvv[conjugate] = nyd - numpy.abs(vvv[conjugate])
    negv = (vvv < 0)
    bugstep = numpy.abs(vvv[negv])
    vvv[negv] = nyd - bugstep

    # numpy bug(?): vectors with values close to zero are handled poorly
    checkhigh = vvv == nyd
    vvv[checkhigh] = nyd - 1
    #print(uuu)
    #print(vvv)

    nvis = vvv.size
    wu = numpy.zeros([width, nvis])
    wv = numpy.zeros([width, nvis])
    step = (ngcf - 1) // width
    rv = vvv - numpy.floor(vvv)
    ru = uuu - numpy.floor(uuu)
    uuu = numpy.floor(uuu).astype(int)
    vvv = numpy.floor(vvv).astype(int)
    p  = ngcf // 2 - numpy.around(step * rv)
    q  = ngcf // 2 - numpy.around(step * ru)
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

    uuun2 = uuu - 2
    uuun1 = uuu - 1
    uuup1 = uuu + 1
    uuup2 = uuu + 2
    uuup3 = uuu + 3
    vvvn2 = vvv - 2
    vvvn1 = vvv - 1
    vvvp1 = vvv + 1
    vvvp2 = vvv + 2
    vvvp3 = vvv + 3

    toohigh = vvvp1 >= nyd
    vvvp1[toohigh] = vvvp1[toohigh] - nyd
    toohigh = vvvp2 >= nyd
    vvvp2[toohigh] = vvvp2[toohigh] - nyd
    toohigh = vvvp3 >= nyd
    vvvp3[toohigh] = vvvp3[toohigh] - nyd
    toolow = uuun2 < 0
    uuun2[toolow] = nyd + (uuun2[toolow])
    toolow = uuun1 < 0
    uuun1[toolow] = nyd + (uuun1[toolow])
    toolow = vvvn2 < 0
    vvvn2[toolow] = nyd + (vvvn2[toolow])
    toolow = vvvn1 < 0
    vvvn1[toolow] = nyd + (vvvn1[toolow])
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

    #print(wu.sum(axis=0).size)
    wu0 = wu[0, :] * Grd[uuun2, vvvn2]
    wu1 = wu[1, :] * Grd[uuun1, vvvn2]
    wu2 = wu[2, :] * Grd[uuu, vvvn2]
    wu3 = wu[3, :] * Grd[uuup1, vvvn2]
    wu4 = wu[4, :] * Grd[uuup2, vvvn2]
    wu5 = wu[5, :] * Grd[uuup3, vvvn2]
    wv0 = wv[0, :] * (wu0 + wu1 + wu2 + wu3 + wu4 + wu5)
    wu0 = wu[0, :] * Grd[uuun2, vvvn1]
    wu1 = wu[1, :] * Grd[uuun1, vvvn1]
    wu2 = wu[2, :] * Grd[uuu, vvvn1]
    wu3 = wu[3, :] * Grd[uuup1, vvvn1]
    wu4 = wu[4, :] * Grd[uuup2, vvvn1]
    wu5 = wu[5, :] * Grd[uuup3, vvvn1]
    wv1 = wv[1, :] * (wu0 + wu1 + wu2 + wu3 + wu4 + wu5)
    wu0 = wu[0, :] * Grd[uuun2, vvv]
    wu1 = wu[1, :] * Grd[uuun1, vvv]
    wu2 = wu[2, :] * Grd[uuu, vvv]
    wu3 = wu[3, :] * Grd[uuup1, vvv]
    wu4 = wu[4, :] * Grd[uuup2, vvv]
    wu5 = wu[5, :] * Grd[uuup3, vvv]
    wv2 = wv[2, :] * (wu0 + wu1 + wu2 + wu3 + wu4 + wu5)
    wu0 = wu[0, :] * Grd[uuun2, vvvp1]
    wu1 = wu[1, :] * Grd[uuun1, vvvp1]
    wu2 = wu[2, :] * Grd[uuu, vvvp1]
    wu3 = wu[3, :] * Grd[uuup1, vvvp1]
    wu4 = wu[4, :] * Grd[uuup2, vvvp1]
    wu5 = wu[5, :] * Grd[uuup3, vvvp1]
    wv3 = wv[3, :] * (wu0 + wu1 + wu2 + wu3 + wu4 + wu5)
    wu0 = wu[0, :] * Grd[uuun2, vvvp2]
    wu1 = wu[1, :] * Grd[uuun1, vvvp2]
    wu2 = wu[2, :] * Grd[uuu, vvvp2]
    wu3 = wu[3, :] * Grd[uuup1, vvvp2]
    wu4 = wu[4, :] * Grd[uuup2, vvvp2]
    wu5 = wu[5, :] * Grd[uuup3, vvvp2]
    wv4 = wv[4, :] * (wu0 + wu1 + wu2 + wu3 + wu4 + wu5)
    wu0 = wu[0, :] * Grd[uuun2, vvvp3]
    wu1 = wu[1, :] * Grd[uuun1, vvvp3]
    wu2 = wu[2, :] * Grd[uuu, vvvp3]
    wu3 = wu[3, :] * Grd[uuup1, vvvp3]
    wu4 = wu[4, :] * Grd[uuup2, vvvp3]
    wu5 = wu[5, :] * Grd[uuup3, vvvp3]
    wv5 = wv[5, :] * (wu0 + wu1 + wu2 + wu3 + wu4 + wu5)
    #print(wu, wv)
    Intp = (wv0 + wv1 + wv2 + wv3 + wv4 + wv5) / w
    #Intp = \
    #    wv[0, :] * ( \
    #    wu[0, :] * Grd[uuu - 2, vvv - 2] + \
    #    wu[1, :] * Grd[uuu - 1, vvv - 2] + \
    #    wu[2, :] * Grd[uuu, vvv - 2] + \
    #    wu[3, :] * Grd[uuu + 1, vvv - 2] + \
    #    wu[4, :] * Grd[uuu + 2, vvv - 2] + \
    #    wu[5, :] * Grd[uuu + 3, vvv - 2]) + \
    #    wv[1, :] * ( \
    #    wu[0, :] * Grd[uuu - 2, vvv - 1] + \
    #    wu[1, :] * Grd[uuu - 1, vvv - 1] + \
    #    wu[2, :] * Grd[uuu, vvv - 1] + \
    #    wu[3, :] * Grd[uuu + 1, vvv - 1] + \
    #    wu[4, :] * Grd[uuu + 2, vvv - 1] + \
    #    wu[5, :] * Grd[uuu + 3, vvv - 1]) + \
    #    wv[2, :] * ( \
    #    wu[0, :] * Grd[uuu - 2, vvv] + \
    #    wu[1, :] * Grd[uuu - 1, vvv] + \
    #    wu[2, :] * Grd[uuu, vvv] + \
    #    wu[3, :] * Grd[uuu + 1, vvv] + \
    #    wu[4, :] * Grd[uuu + 2, vvv] + \
    #    wu[5, :] * Grd[uuu + 3, vvv]) + \
    #    wv[3, :] * ( \
    #    wu[0, :] * Grd[uuu - 2, vvv + 1] + \
    #    wu[1, :] * Grd[uuu - 1, vvv + 1] + \
    #    wu[2, :] * Grd[uuu, vvv + 1] + \
    #    wu[3, :] * Grd[uuu + 1, vvv + 1] + \
    #    wu[4, :] * Grd[uuu + 2, vvv + 1] + \
    #    wu[5, :] * Grd[uuu + 3, vvv + 1]) + \
    #    wv[4, :] * ( \
    #    wu[0, :] * Grd[uuu - 2, vvv + 2] + \
    #    wu[1, :] * Grd[uuu - 1, vvv + 2] + \
    #    wu[2, :] * Grd[uuu, vvv + 2] + \
    #    wu[3, :] * Grd[uuu + 1, vvv + 2] + \
    #    wu[4, :] * Grd[uuu + 2, vvv + 2] + \
    #    wu[5, :] * Grd[uuu + 3, vvv + 2]) + \
    #    wv[5, :] * ( \
    #    wu[0, :] * Grd[uuu - 2, vvv + 3] + \
    #    wu[1, :] * Grd[uuu - 1, vvv + 3] + \
    #    wu[2, :] * Grd[uuu, vvv + 3] + \
    #    wu[3, :] * Grd[uuu + 1, vvv + 3] + \
    #    wu[4, :] * Grd[uuu + 2, vvv + 3] + \
    #    wu[5, :] * Grd[uuu + 3, vvv + 3])

    #for i in xrange(nvis):

    #    u = uuu[i]
    #    v = vvv[i]
        #conjugate = (uuu < 0)

        #if conjugate:
    #    if u < 0:
    #      u = -1 * uuu[i]
    #      v = -1 * vvv[i]
        #else:
          #u = uuu
          #v = vvv

    #    if v < 0:
    #        v = nyd - numpy.abs(v)

    #    ju = numpy.int(u)
    #    jv = numpy.int(v)
        #p  = ngcf / 2 #- numpy.around(step * (v - jv))
        #q  = ngcf / 2 #- numpy.around(step * (u - ju))
        #print(ju, jv, p, q, step)
        #print(numpy.around(step * (u - ju)), numpy.around(step * (u - numpy.int(u))))

        # - from here width = 6 (see ModGrid in Model.for)
        #wu1 = gcf[q - 2 * step]
        #wu2 = gcf[q - step]
        #wu3 = gcf[q]
        #wu4 = gcf[q + step]
        #wu5 = gcf[q + 2 * step]
        ##wu6 = gcf[q + 3 * step]
        #wuindx = range(q - 2 * step, q + 3 * step, step)
        #wu = gcf[wuindx]
        #wvindx = range(p - 2 * step, p + 3 * step, step)
        #wv = gcf[wvindx]

        #wv1 = gcf[p - 2 * step]
        #wv2 = gcf[p - step]
        #wv3 = gcf[p]
        #wv4 = gcf[p + step]
        #wv5 = gcf[p + 2 * step]
        ##wv6 = gcf[p + 3 * step]

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

    #    jup1 = ju + 1
    #    jup2 = ju + 2
        #jup3 = ju + 3
    #    if jv + 1 >= nyd:
    #        jvp1 = jv + 1 - nyd
    #    else:
    #        jvp1 = jv + 1 
    #    if jv+2 >= nyd:
    #        jvp2 = jv + 2 - nyd
    #    else:
    #        jvp2 = jv + 2 
        #if jv + 3 >= nyd:
        #    jvp3 = jv + 3 - nyd
        #else:
        #    jvp3 = jv + 3

        ##w = (wu1 + wu2 + wu3 + wu4 + wu5 + wu6) * \
        ##        (wv1 + wv2 + wv3 + wv4 + wv5 + wv6)
        #w = (wu1 + wu2 + wu3 + wu4 + wu5) * \
        #        (wv1 + wv2 + wv3 + wv4 + wv5)

        #Intpi = \
        #    wv1 * ( wu1 * Grd[jun2, jvn2] + wu2 * Grd[jun1, jvn2] + \
        #       wu3 * Grd[ju, jvn2] + wu4 * Grd[jup1, jvn2] + \
        #       wu5 * Grd[jup2, jvn2] + wu6 * Grd[jup3, jvn2]) + \
        #    wv2 * (wu1 * Grd[jun2, jvn1] + wu2 * Grd[jun1, jvn1] + \
        #       wu3 * Grd[ju, jvn1] + wu4 * Grd[jup1, jvn1] + \
        #       wu5 * Grd[jup2, jvn1] + wu6 * Grd[jup3, jvn1]) + \
        #    wv3 * (wu1 * Grd[jun2, jv] + wu2 * Grd[jun1, jv] + \
        #       wu3 * Grd[ju, jv] + wu4 * Grd[jup1, jv] + \
        #       wu5 * Grd[jup2, jv] + wu6 * Grd[jup3, jv]) + \
        #    wv4 * (wu1 * Grd[jun2, jvp1] + wu2 * Grd[jun1, jvp1] + \
        #       wu3 * Grd[ju, jvp1] + wu4 * Grd[jup1, jvp1] + \
        #       wu5 * Grd[jup2, jvp1] + wu6 * Grd[jup3, jvp1]) + \
        #    wv5 * (wu1 * Grd[jun2, jvp2] + wu2 * Grd[jun1, jvp2] + \
        #       wu3 * Grd[ju, jvp2] + wu4 * Grd[jup1, jvp2] + \
        #       wu5 * Grd[jup2, jvp2] + wu6 * Grd[jup3, jvp2]) + \
        #    wv6 * (wu1 * Grd[jun2, jvp3] + wu2 * Grd[jun1, jvp3] + \
        #       wu3 * Grd[ju, jvp3] + wu4 * Grd[jup1, jvp3] + \
        #       wu5 * Grd[jup2, jvp3] + wu6 * Grd[jup3, jvp3])
    #    Intpi = \
    #        wv[0] * ( wu[0] * Grd[jun2, jvn2] + wu[1] * Grd[jun1, jvn2] + \
    #           wu[2] * Grd[ju, jvn2] + wu[3] * Grd[jup1, jvn2] + \
    #           wu[4] * Grd[jup2, jvn2]) + \
    #        wv[1] * (wu[0] * Grd[jun2, jvn1] + wu[1] * Grd[jun1, jvn1] + \
    #           wu[2] * Grd[ju, jvn1] + wu[3] * Grd[jup1, jvn1] + \
    #           wu[4] * Grd[jup2, jvn1]) + \
    #        wv[2] * (wu[0] * Grd[jun2, jv] + wu[1] * Grd[jun1, jv] + \
    #           wu[2] * Grd[ju, jv] + wu[3] * Grd[jup1, jv] + \
    #           wu[4] * Grd[jup2, jv]) + \
    #        wv[3] * (wu[0] * Grd[jun2, jvp1] + wu[1] * Grd[jun1, jvp1] + \
    #           wu[2] * Grd[ju, jvp1] + wu[3] * Grd[jup1, jvp1] + \
    #           wu[4] * Grd[jup2, jvp1]) + \
    #        wv[4] * (wu[0] * Grd[jun2, jvp2] + wu[1] * Grd[jun1, jvp2] + \
    #           wu[2] * Grd[ju, jvp2] + wu[3] * Grd[jup1, jvp2] + \
    #           wu[4] * Grd[jup2, jvp2])

    #    print(u, ju, v, jv)
        #cdef double wu1G1 = wu1 * Grd[jun2, jvn2]
        #cdef double wu2G1 = wu2 * Grd[jun1, jvn2]
        #cdef double wu3G1 = wu3 * Grd[ju, jvn2]
        #cdef double wu4G1 = wu4 * Grd[jup1, jvn2]
        #cdef double wu5G1 = wu5 * Grd[jup2, jvn2]
        #cdef double wu6G1 = wu6 * Grd[jup3, jvn2]
        #cdef double wv1wu = wv1 * (wu1G1 + wu2G1 + wu3G1 + wu4G1 + wu5G1 + wu6G1)
        #cdef double wu1G2 = wu1 * Grd[jun2, jvn1]
        #cdef double wu2G2 = wu2 * Grd[jun1, jvn1]
        #cdef double wu3G2 = wu3 * Grd[ju, jvn1]
        #cdef double wu4G2 = wu4 * Grd[jup1, jvn1]
        #cdef double wu5G2 = wu5 * Grd[jup2, jvn1]
        #cdef double wu6G2 = wu6 * Grd[jup3, jvn1]
        #cdef double wv2wu = wv2 * (wu1G2 + wu2G2 + wu3G2 + wu4G2 + wu5G2 + wu6G2)
        #cdef double wu1G3 = wu1 * Grd[jun2, jv]
        #cdef double wu2G3 = wu2 * Grd[jun1, jv]
        #cdef double wu3G3 = wu3 * Grd[ju, jv]
        #cdef double wu4G3 = wu4 * Grd[jup1, jv]
        #cdef double wu5G3 = wu5 * Grd[jup2, jv]
        #cdef double wu6G3 = wu6 * Grd[jup3, jv]
        #cdef double wv3wu = wv3 * (wu1G3 + wu2G3 + wu3G3 + wu4G3 + wu5G3 + wu6G3)
        #cdef double wu1G4 = wu1 * Grd[jun2, jvp1]
        #cdef double wu2G4 = wu2 * Grd[jun1, jvp1]
        #cdef double wu3G4 = wu3 * Grd[ju, jvp1]
        #cdef double wu4G4 = wu4 * Grd[jup1, jvp1]
        #cdef double wu5G4 = wu5 * Grd[jup2, jvp1]
        #cdef double wu6G4 = wu6 * Grd[jup3, jvp1]
        #cdef double wv4wu = wv4 * (wu1G4 + wu2G4 + wu3G4 + wu4G4 + wu5G4 + wu6G4)
        #cdef double wu1G5 = wu1 * Grd[jun2, jvp2]
        #cdef double wu2G5 = wu2 * Grd[jun1, jvp2]
        #cdef double wu3G5 = wu3 * Grd[ju, jvp2]
        #cdef double wu4G5 = wu4 * Grd[jup1, jvp2]
        #cdef double wu5G5 = wu5 * Grd[jup2, jvp2]
        #cdef double wu6G5 = wu6 * Grd[jup3, jvp2]
        #cdef double wv5wu = wv5 * (wu1G5 + wu2G5 + wu3G5 + wu4G5 + wu5G5 + wu6G5)
        #cdef double wu1G6 = wu1 * Grd[jun2, jvp3]
        #cdef double wu2G6 = wu2 * Grd[jun1, jvp3]
        #cdef double wu3G6 = wu3 * Grd[ju, jvp3]
        #cdef double wu4G6 = wu4 * Grd[jup1, jvp3]
        #cdef double wu5G6 = wu5 * Grd[jup2, jvp3]
        #cdef double wu6G6 = wu6 * Grd[jup3, jvp3]
        #cdef double wv6wu = wv6 * (wu1G6 + wu2G6 + wu3G6 + wu4G6 + wu5G6 + wu6G6)

        #cdef double Intp = wv1wu + wv2wu + wv3wu + wv4wu + wv5wu + wv6wu

    #  Conjugate the data if necessary.
    Intp[conjugate] = numpy.conjugate(Intp[conjugate])

    #    Intp[i] = Intpi / w

    #print(Intp)
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
    ln2 = numpy.log(2)
    nxd = 2 ** numpy.ceil(numpy.log(2. * nx) / ln2)	# - padding w/ zeros
    nyd = 2 ** numpy.ceil(numpy.log(2. * ny) / ln2)	# - padding w/ zeros
    nxd = numpy.long(nxd)
    nyd = numpy.long(nyd)
    dx = modelheader['CDELT1'] * pi / 180  # - pixel size in radians
    dy = modelheader['CDELT2'] * pi / 180  # - pixel size in radians 
    du = 1. / (dx * nxd)		# - size of u cell
    dv = 1. / (dy * nyd)		# - size of v cell
    nu = nxd			# - # of u cells
    nv = nyd			# - # of v cells
    u0 = nu / 2 + 1
    v0 = nv / 2 + 1
    umax = 0.5 * (nxd - 1 - width)
    vmax = 0.5 * (nyd - 1 - width)

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
        numpy.cos(pcd[1] * pi / 180.) * pi / 180.
    decref2 = (pcd[1] - pcm[1]) * pi / 180.

    # - CALCULATE UVW CONVERSION MATRIX FOR GEOMETRY DIFFERENCES
    # phase center for the model
    ucoeff, vcoeff, wcoeff = grid.coGeom(pcm, pcd)

    # - shift array and pad with zeros
    shifted = numpy.zeros([nyd, nxd])
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
    image = numpy.roll(image, ny // 2, axis=0)
    image = numpy.roll(image, nx // 2, axis=1)

    #mu_grid, dump = numpy.meshgrid(mu, numpy.zeros(nxd) + 1)
    #mu_grid_shifty = numpy.roll(mu_grid, -1 * u0int, axis=0)
    #mu_grid_shiftyx = numpy.roll(mu_grid, -1 * u0int, axis=1)
    #foo = du * mu_grid_shiftyx
    #tmu = numpy.transpose(mu)
    #tmu_grid, dump = numpy.meshgrid(tmu, numpy.zeros(nxd) + 1)
    #tmu_grid_shifty = numpy.roll(tmu_grid, -1 * u0int, axis=0)
    #tmu_grid_shiftyx = numpy.roll(tmu_grid, -1 * u0int, axis=1)
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
        mcorr = numpy.outer(ycorr, xcorr)
        image = image / mcorr
        #time_modgrid = time.time()-start
        #print('time to run ModCorr: ', time_modgrid, 'seconds')

    # - take FFT
    #print('ModPlane (model.for)')
    #cdef numpy.ndarray mvis_real
    #cdef numpy.ndarray mvis_imag
    mvis = numpy.fft.fft2(image)
    mvis = numpy.conjugate(mvis)

    # - Follow ModGrid in model.for
    #print('ModGrid (model.for)')
    #cdef numpy.ndarray gcf
    ngcf = width * ((maxgcf - 1) // width) + 1
    gcf = grid.gcffun(ngcf, width, alpha)

    #cdef numpy.ndarray uu
    #cdef numpy.ndarray vv

    # - interpolate grid
    if igrid == 1:
        #print('interpolating grid')
        # mvis = transpose(mvis)
        uu = ud / du  # calculate fractional index in the grid
        vv = vd / dv  # calculate fractional index in the grid
        #mvis_opt = numpy.zeros(nvis) + 1.j * numpy.zeros(nvis)
        # - weighted interpolation
        #import pyximport
        #pyximport.install()
        #import ModGrid
        #for i in numpy.arange(nvis):
        #    if (numpy.abs(uu[i]) > umax) or (numpy.abs(vv[i]) > vmax):
        #        mvis_opt[i] = 0.
        #    else:
        #        start = time.time()
        #mvis_real = numpy.real(mvis)
        #mvis_imag = numpy.imag(mvis)
        #uuint = numpy.around(uu).astype(int)
        #vvint = numpy.around(vv).astype(int)
        #ru = uu - uuint
        #rv = vv - vvint
        mvis_opt = ModGrid1(uu, vv, u0, v0, mvis, gcf, ngcf, nyd, width)
        #print(uu)
        #print(vv)
        overmax = numpy.abs(vv) > vmax
        mvis_opt[overmax] = 0.
        overmax = numpy.abs(uu) > umax
        mvis_opt[overmax] = 0.
        #mvis_opt1 = numpy.zeros(nvis) + 1.j * numpy.zeros(nvis)
        #p = numpy.zeros(nvis)
        #q = numpy.zeros(nvis)
        #for i in numpy.arange(nvis):
        #    if (numpy.abs(uu[i]) > umax) or (numpy.abs(vv[i]) > vmax):
        #        mvis_opt1[i] = 0.
        #    else:
        #        start = time.time()
        #        mvis_real = numpy.real(mvis)
        #        mvis_imag = numpy.imag(mvis)
        #        mvis_opt1[i], p[i], q[i] = ModGrid.ModGrid(uu[i], vv[i], mvis, 1, 1, 1, \
        #                gcf, ngcf, nyd)
        #        time_modgrid = time_modgrid + time.time()-start
        #        start = time.time()
                #mvis_opt[i] = grid.ModShift(ud[i], vd[i], raref1, decref1, raref2, \
                #        decref2, 1, 1, mvis_opt[i])

        #print(uu)
        #diffmvis = numpy.abs((mvis_opt - mvis_opt1) / mvis_opt)
        #print(diffmvis)
        #print(diffmvis.min(), diffmvis.max())
        #mvis_opt_imag = ModGrid(uuint, vvint, mvis_imag, gcf, ngcf, nyd, width, step)
        #mvis_opt = mvis_opt_real + 1.j * mvis_opt_imag
        #        time_modgrid = time_modgrid + time.time()-start
        #        start = time.time()
        mvis_opt = grid.ModShift(ud, vd, raref1, decref1, raref2, \
                        decref2, 1, 1, mvis_opt)
        #for i in numpy.arange(nvis):
        #    if (numpy.abs(uu[i]) > umax) or (numpy.abs(vv[i]) > vmax):
        #        mvis_opt[i] = 0.
        #    else:
        #        mvis_opt[i] = grid.ModShift(ud[i], vd[i], raref1, decref1, raref2, \
        #                decref2, 1, 1, mvis_opt[i])
                #time_modshift = time_modshift + time.time()-start
        #print(time_modgrid, time_modshift)
    if igrid == 2:
        u0 = nxd // 2 + 1		# - central u cell id
        #v0 = nyd // 2 + 1		# - central v cell id
        mu = numpy.arange(nxd)
        u0int = numpy.int(u0)
        mu[u0int:] = u0int - nxd + numpy.arange(u0int - 2)
        mu1d = du * mu
        import time
        # - Nearest neighbor search
        uu = ud / du #  calculate index in the grid
        vv = vd / dv #  calculate index in the grid
        mvis_opt = numpy.zeros(nvis) + 1.j * numpy.zeros(nvis)
        mv1d = -1 * mu1d
        imu1d = numpy.argsort(mu1d)
        smu1d = mu1d[imu1d]
        imv1d = numpy.argsort(mv1d)
        smv1d = mv1d[imv1d]
        time_modgrid = 0
        for i in numpy.arange(nvis):
            start = time.time()
            iu = (numpy.abs(smu1d - u[i])).argmin()
            iv = (numpy.abs(smv1d - v[i])).argmin()
            mvis_opt[i] = mvis[imv1d[iv], imu1d[iu]]
            mvis_opt[i] = grid.ModShift(ud[i], vd[i], raref1, decref1, raref2, \
                    decref2, 1, 1, mvis_opt[i])
            time_modgrid = time_modgrid + time.time()-start
        #print(time_modgrid)

    #mvis_opt = numpy.conjugate(mvis_opt)
    # - read in visibilities MIRIAD produces
    #imvis = pyfits.getdata(imvfile)
    #mirvisheader = pyfits.getheader(imvfile)
    #ruv = numpy.sqrt(u ** 2 + v ** 2)

    # read c version
    #readcol,'../myvis.dat',cu,cv,cr,ci

    # - plot the results
    #rvis = imvis['DATA'][:, 0, 0, 0, 0, 0]
    #ivis = imvis['DATA'][:, 0, 0, 0, 0, 1]
    #wgt = imvis['DATA'][:, 0, 0, 0, 0, 2]
    #ampvis = numpy.sqrt(rvis ** 2 + ivis ** 2)

    # Real part
    #mvis_amp = numpy.sqrt(numpy.real(mvis_opt) ** 2 + numpy.imag(mvis_opt) ** 2)
    #mvis_amp = numpy.abs(mvis_opt)
    #plt.clf()
    #plt.plot(ruv, mvis_amp-ampvis, ',', color='red', zorder=2)
    #plt.axis([0, 300, 0, 1])
    #plt.plot(ruv, numpy.real(mvis_opt), 'o', ms=2, mec='red', mfc='red')  # Python
    #plt.plot(ruv, rvis, '.', ms=1, color='black')  # MIRIAD
    #plt.show()
    #saveloc = 'real_python_miriad'
    #savefig(saveloc)

     # Imaginary part
    #import matplotlib.pyplot as plt
    #from pylab import savefig
    #plt.clf()
    #plt.plot(ruv, numpy.abs(mvis_opt), 'o', ms=2, mec='red', mfc='red')  # Python
    #plt.loglog()
    ##plt.plot(ruv, ivis, 'o', mfc='black', ms=1)  # MIRIAD
    #saveloc = 'imaginary_python_miriad'
    #savefig(saveloc)

    # calculate chi-squared 
    #chi2v = wgt * (numpy.real(mvis_opt) - rvis) ** 2 + \
    #    wgt * (numpy.imag(mvis_opt) - ivis) ** 2
    #print('chi squared value is :', chi2v.sum())
    # 4.7212017638007328e-06
    #pdb.set_trace()
    #print(numpy.abs(mvis_opt).max())
    #print(numpy.abs(mvis_opt).min())
    #import pdb; pdb.set_trace()
    return mvis_opt
