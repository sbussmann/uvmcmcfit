# Plot results from emcee
# Shane Bussmann
# 2012 March 31

def plotmodel(smacut, modcut, rms, modradx, modrady, celldata, outcrit):

    import math
    import numpy

    diffcut = smacut - modcut

    #--------------------------------------------------------------------------

    # Step 8: Draw the figures
    #import aplpy
    import matplotlib.pyplot as mpl


    fig = mpl.figure(figsize=(11, 5))
    ax = fig.add_subplot(1, 2, 1, aspect='equal')
    mpl.subplots_adjust(left=0.10, right=0.98, top=0.99, bottom=0.10, wspace=0.3)

    inv_mod = modcut.max()- modcut
    modx = ( numpy.arange(2*modradx+2) - modradx ) * celldata
    mody = ( numpy.arange(2*modrady+2) - modrady ) * celldata
    mpl.pcolor(modx, mody, inv_mod, cmap='gray', edgecolor='None')
    mpl.pcolor(modx, mody, inv_mod, cmap='gray', edgecolor='none')
    mpl.pcolor(modx, mody, inv_mod, cmap='gray', edgecolor='none')
    #mpl.imshow(inv_mod, cmap='gray')

    levs = numpy.array([-4,-2,2,4,6,8,10,12,14,16,18,20,22])*rms
    cline = ['dashed','dashed','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid']
    cmodx = ( numpy.arange(2*modradx+1) - modradx ) * celldata
    cmody = ( numpy.arange(2*modrady+1) - modrady ) * celldata
    mpl.contour(cmodx, cmody, smacut, colors='red', levels=levs, linestyles=cline, \
        linewidths=1.5)

    mpl.minorticks_on()
    mpl.tick_params(width=2, which='both')
    mpl.tick_params(length=2, which='minor')
    mpl.tick_params(length=4, which='major')
    mpl.xlabel(r'$\Delta$RA (arcsec)')
    mpl.ylabel(r'$\Delta$Dec (arcsec)')

    # read in critical curve data
    crit_data = asciitable.read(outcrit, Reader=asciitable.NoHeader)

    cx1 = crit_data['col1']
    cy1 = crit_data['col2']
    cu1 = crit_data['col3']
    cv1 = crit_data['col4']
    cx2 = crit_data['col5']
    cy2 = crit_data['col6']
    cu2 = crit_data['col7']
    cv2 = crit_data['col8']

    # overplot caustics and critical curves
    nclines = cx1.size
    #for i=0,nclines-85 do begin
    for i in numpy.arange(0,nclines):
        xarr = numpy.array([cu1[i], cu2[i]]) #- offx*celldata
        yarr = numpy.array([cv1[i], cv2[i]]) #- offy*celldata
        mpl.plot(xarr, yarr, color='cyan', lw=2)
        xarr = numpy.array([cx1[i], cx2[i]]) #- offx*celldata
        yarr = numpy.array([cy1[i], cy2[i]]) #- offy*celldata
        mpl.plot(xarr, yarr, color='orange', lw=2)

    # overplot position of the background source and the foreground lenses
    xxx = numpy.array( bestfit[0][3])#-offx*celldata )
    yyy = numpy.array( bestfit[0][4])#-offy*celldata )
    mpl.plot(xxx, yyy, 'o', color='blue', fillstyle='full', \
        mec='white', mew=0., ms=9., label='Source Position')
    xxx1 = numpy.array( 0)#-offx*celldata )
    yyy1 = numpy.array( 0)#-offy*celldata )
    mpl.plot(xxx1, yyy1, '+', ms=9., mfc='black', mec='black', mew=2., \
        label='Lens Position')

    import matplotlib.font_manager as fm
    prop = fm.FontProperties(size='medium')
    mpl.legend(numpoints=1, prop=prop, handlelength=1.0)
    #mpl.legend(numpoints=1, handletextpad=0.05, labelspacing=0.25, handlelength=1.5, \
    #   borderpad=0.25, markerscale=1.0)

    #ax = fig.add_subplot(2,2,i+1)
    ax.text(0.05, 0.95, 'G9-19', transform=ax.transAxes, \
        fontsize='large', va='top')
    ax.text(0.05, 0.85, r'$n_s =  0.5$', transform=ax.transAxes, \
        fontsize='large', va='top')
    #ax.text(0.05, 0.20, r'$\mu = 3.56$', transform=ax.transAxes, \
    #    fontsize='large', va='top')
    #ax.text(0.05, 0.13, r'$\chi_\nu^2 =  0.90$', transform=ax.transAxes, \
    #    fontsize='large', va='top')

    #from matplotlib.patches import Ellipse
    #print 
    #e = Ellipse((0.88,0.10), bmin/5., bmaj/5., angle=bpa, ec='black', hatch='/', \
    #    lw=2, transform=ax.transAxes, fc='white')
    #ax.add_artist(e)

    #mpl.text(0.1, 0.8, 'G15-141')

    # set some legend properties.  All the code below is optional.  The
    # defaults are usually sensible but if you need more control, this
    # shows you how
    #leg = mpl.gca().get_legend()
    #ltext  = leg.get_texts()  # all the text.Text instance in the legend
    #llines = leg.get_lines()  # all the lines.Line2D instance in the legend
    #frame  = leg.get_frame()  # the patch.Rectangle instance surrounding the legend

    # see text.Text, lines.Line2D, and patches.Rectangle for more info on
    # the settable properties of lines, text, and rectangles
    #frame.set_facecolor('0.80')      # set the frame face color to light gray
    #mpl.setp(ltext, fontsize='medium')    # the legend text fontsize
    #mpl.setp(llines, linewidth=1.0)      # the legend linewidth
    #leg.draw_frame(False)           # don't draw the legend frame

    #mpl.text(0.1, 0.8, objname, transform = ax.transAxes)

    mpl.subplot(1, 2, 2, aspect='equal')

    inv_mod = diffcut.max()- diffcut
    modx = ( numpy.arange(2*modradx+2) - modradx ) * celldata
    mody = ( numpy.arange(2*modrady+2) - modrady ) * celldata
    mpl.pcolor(modx, mody, inv_mod, cmap='gray', edgecolor='None')
    mpl.pcolor(modx, mody, inv_mod, cmap='gray', edgecolor='none')
    mpl.pcolor(modx, mody, inv_mod, cmap='gray', edgecolor='none')
    #mpl.imshow(inv_mod, cmap='gray')

    levs = numpy.array([-4,-3,-2,-1,1,2,3,4])*rms
    cline = ['dashed','dashed','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid']
    cmodx = ( numpy.arange(2*modradx+1) - modradx ) * celldata
    cmody = ( numpy.arange(2*modrady+1) - modrady ) * celldata
    mpl.contour(cmodx, cmody, diffcut, colors='black', levels=levs, linewidths=1.5)

    mpl.minorticks_on()
    mpl.tick_params(width=2, which='both')
    mpl.tick_params(length=2, which='minor')
    mpl.tick_params(length=4, which='major')
    mpl.xlabel(r'$\Delta$RA (arcsec)')
    mpl.ylabel(r'$\Delta$Dec (arcsec)')

    ndof = (2*modradx+1) * (2*modrady+1)
    #print numpy.sum(diffcut**2)/rms**2, ndof, numpy.sum(diffcut**2)/rms**2/ndof

    #f1 = aplpy.FITSFigure(modcutloc, figure=fig, subplot=[0.1,0.1,0.35,0.8])
    ##f1.set_tick_labels_font(size='x-small')
    ##f1.set_axis_labels_font(size='small')
    #f1.show_grayscale(invert=True, vmin=0)
    #axes(subplot=[0.1,0.1,0.35,0.8])
    #
    ##f1a = aplpy.FITSFigure(smacutloc, figure=fig, subplot=[0.1,0.1,0.35,0.8])
    #f1.show_contour(smacutloc, levels=levs, colors='red', linestyles=cline, \
        #linewidths=2)
    #
    #f1.add_beam()
    #f1.beam.show()
    #f1.beam.set(edgecolor='black', lw=3, hatch='/', facecolor='white')
    #
    #f2 = aplpy.FITSFigure(diffcutloc, figure=fig, subplot=[0.5,0.1,0.35,0.8])
    ##f2.set_tick_labels_font(size='x-small')
    ##f2.set_axis_labels_font(size='small')
    #f2.show_grayscale()
    #
    ##f2.hide_yaxis_label()
    ##f2.hide_ytick_labels()
    #
    #fig.canvas.draw()

    from pylab import savefig
    savefig('tex_demo')
