"""

2014 July 10

Shane Bussmann

Convert config.py to config.yaml

"""

import config
import yaml


def grabvalues(config, parname):
    priorparname = 'Prior_' + parname
    prior = getattr(config, priorparname)
    initparname = 'Init_' + parname
    init = getattr(config, initparname)
    limits = [prior[0], init[0], init[1], prior[1]]
    fixed = prior[2]
    print(fixed)
    if fixed != 'free':
        fixed = fixed.split("_")
        fixed.reverse()
        fixed = " ".join(fixed)
    priorshape = prior[3]
    return limits, fixed, priorshape

objectname = config.ObjectName
ImageName = config.ImageName
opticalimage = config.OpticalImage
opticaltag = config.OpticalTag
Nwalkers = config.Nwalkers
Nthreads = config.Nthreads
lnLike = config.lnLike
uvdata = config.VisFile
parallel = config.ParallelProcessingMode
if parallel == 'MPI':
    MPI = True
else:
    MPI = False

configdictlist = [
        ('ImageName', ImageName),
        ('LogLike', lnLike),
        ('MPI', MPI),
        ('Nthreads', Nthreads),
        ('Nwalkers', Nwalkers),
        ('ObjectName', objectname),
        ('OpticalImage', opticalimage),
        ('OpticalTag', opticaltag),
        ('UVData', uvdata)
        ]

regionid = config.RegionID
nregions = len(regionid)

for iregion in range(nregions):
    sr = str(iregion)
    RACentroid = config.RACentroid[iregion]
    DecCentroid = config.DecCentroid[iregion]
    Oversample = config.Oversample[iregion]
    RadialExtent = config.RadialExtent[iregion]
    regionkey = 'Region' + sr
    regiondictlist = [
            ('RACentroid', RACentroid),
            ('DecCentroid', DecCentroid),
            ('Oversample', Oversample),
            ('RadialExtent', RadialExtent)
            ]
    Nlens = config.Nlens[iregion]
    for ilens in range(Nlens):
        sl = str(ilens)
        lenskey = 'Lens' + sl

        # Einstein Radius
        parname = 'EinsteinRadius_Lens' + sl + '_Region' + sr
        limits, fixed, priorshape = grabvalues(config, parname)
        EinsteinRadius = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # DeltaRA
        parname = 'DeltaRA_Lens' + sl + '_Region' + sr
        limits, fixed, priorshape = grabvalues(config, parname)
        DeltaRA = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # DeltaDec
        parname = 'DeltaDec_Lens' + sl + '_Region' + sr
        limits, fixed, priorshape = grabvalues(config, parname)
        DeltaDec = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # AxialRatio
        parname = 'AxialRatio_Lens' + sl + '_Region' + sr
        limits, fixed, priorshape = grabvalues(config, parname)
        AxialRatio = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # PositionAngle
        parname = 'PositionAngle_Lens' + sl + '_Region' + sr
        limits, fixed, priorshape = grabvalues(config, parname)
        PositionAngle = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        lensdictlist = [
                ('EinsteinRadius', EinsteinRadius),
                ('DeltaRA', DeltaRA),
                ('DeltaDec', DeltaDec),
                ('AxialRatio', AxialRatio),
                ('PositionAngle', PositionAngle)
                ]
        lensdict = dict(lensdictlist)
        lenstuple = (lenskey, lensdict)
        regiondictlist.append(lenstuple)

    Nsource = config.Nsource[iregion]
    for isource in range(Nsource):
        ss = str(isource)
        sourcekey = 'Source' + ss

        # Model type
        parname = 'ModelMorphology_Source' + ss + '_Region' + sr
        modeltype = getattr(config, parname)
        if modeltype == 'gaussian':
            modeltype = 'Gaussian'

        # IntrinsicFlux
        parname = 'IntrinsicFlux_Source' + ss + '_Region' + sr
        limits, fixed, priorshape = grabvalues(config, parname)
        IntrinsicFlux = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # DeltaRA
        parname = 'DeltaRA_Source' + ss + '_Region' + sr
        limits, fixed, priorshape = grabvalues(config, parname)
        DeltaRA = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # DeltaDec
        parname = 'DeltaDec_Source' + ss + '_Region' + sr
        limits, fixed, priorshape = grabvalues(config, parname)
        DeltaDec = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # Size
        parname = 'Size_Source' + ss + '_Region' + sr
        limits, fixed, priorshape = grabvalues(config, parname)
        Size = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # AxialRatio
        parname = 'AxialRatio_Source' + ss + '_Region' + sr
        limits, fixed, priorshape = grabvalues(config, parname)
        AxialRatio = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # PositionAngle
        parname = 'PositionAngle_Source' + ss + '_Region' + sr
        limits, fixed, priorshape = grabvalues(config, parname)
        PositionAngle = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        sourcedictlist = [
                ('LightProfile', modeltype),
                ('IntrinsicFlux', IntrinsicFlux),
                ('DeltaRA', DeltaRA),
                ('DeltaDec', DeltaDec),
                ('EffectiveRadius', Size),
                ('AxialRatio', AxialRatio),
                ('PositionAngle', PositionAngle)
                ]
        sourcedict = dict(sourcedictlist)
        sourcetuple = (sourcekey, sourcedict)
        regiondictlist.append(sourcetuple)

    regiondict = dict(regiondictlist)
    regiontuple = (regionkey, regiondict)
    configdictlist.append(regiontuple)

configdict = dict(configdictlist)

configfile = open('config.yaml', 'w')
yaml.dump(configdict, configfile)
