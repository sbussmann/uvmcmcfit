"""

A collection of utility routines to help with running miriad.

"""

def makeScript(parameters):

    """ Make a c-shell script to do imaging. """

    f = open('image.csh', 'w')
    f.write('#! /bin/csh\n')
    f.write('set VIS=' + parameters[0] + '\n')
    f.write('set OUT=' + parameters[1] + '\n')
    f.write('set IMSIZE=' + parameters[2] + '\n')
    f.write('set CELL=' + parameters[3] + '\n')
    f.write('set NITERS=' + parameters[4] + '\n')
    f.write('set CUTOFF=' + parameters[5] + '\n')
    f.write('set CUTOFF2=' + parameters[6] + '\n')
    f.write('set ROBUST=' + parameters[7] + '\n')
    f.write('set REGION=' + parameters[8] + '\n')
    f.write('set REGION2=' + parameters[9] + '\n')
    f.write('set GAIN=' + parameters[10] + '\n')
    f.write('set FWHM=' + parameters[11] + '\n')
    f.write('set SUP=' + parameters[12] + '\n')
    f.write('rm -rf $OUT.map $OUT.beam \n')
    f.write('invert vis=$VIS map=$OUT.map beam=$OUT.beam imsize=$IMSIZE cell=$CELL fwhm=$FWHM sup=$SUP options=systemp,mfs robust=$ROBUST \n')
    f.write('rm -rf $OUT.cc $OUT.cc1 $OUT.cc2 $OUT.cm $OUT.cm1 $OUT.cm2 \n')
    f.write('clean map=$OUT.map beam=$OUT.beam out=$OUT.cc gain=$GAIN niters=$NITERS region=$REGION cutoff=$CUTOFF \n')
    f.write('clean map=$OUT.map beam=$OUT.beam out=$OUT.cc2 gain=$GAIN niters=$NITERS region=$REGION2 cutoff=$CUTOFF2 \n')
    f.write('clean map=$OUT.map beam=$OUT.beam out=$OUT.cc1 gain=$GAIN niters=$NITERS region=$REGION2 cutoff=$CUTOFF2 model=$OUT.cc \n')
    f.write('restor map=$OUT.map model=$OUT.cc beam=$OUT.beam out=$OUT.cm \n')
    f.write('restor map=$OUT.map model=$OUT.cc2 beam=$OUT.beam out=$OUT.cm2 \n')
    f.write('restor map=$OUT.map model=$OUT.cc1 beam=$OUT.beam out=$OUT.cm1 \n')
    f.write('rm -rf $OUT.cm.fits $OUT.fits $OUT.cm2.fits \n')
    f.write('fits in=$OUT.cm1 out=$OUT.fits op=xyout \n')
    f.close()
