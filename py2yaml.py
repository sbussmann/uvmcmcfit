"""

2014 July 10

Shane Bussmann

Convert config.py to config.yaml

"""

import config
import yaml


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
        parname = 'Prior_EinsteinRadius_Lens' + sl + '_Region' + sr
        prior = getattr(config, parname)
        parname = 'Init_EinsteinRadius_Lens' + sl + '_Region' + sr
        init = getattr(config, parname)
        limits = [prior[0], init[0], init[1], prior[1]]
        fixed = prior[2]
        priorshape = prior[3]
        EinsteinRadius = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # DeltaRA
        DeltaRA = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])
        parname = 'Prior_DeltaRA_Lens' + sl + '_Region' + sr
        prior = getattr(config, parname)
        parname = 'Init_DeltaRA_Lens' + sl + '_Region' + sr
        init = getattr(config, parname)
        limits = [prior[0], init[0], init[1], prior[1]]
        fixed = prior[2]
        priorshape = prior[3]
        DeltaRA = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # DeltaDec
        DeltaDec = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])
        parname = 'Prior_DeltaDec_Lens' + sl + '_Region' + sr
        prior = getattr(config, parname)
        parname = 'Init_DeltaDec_Lens' + sl + '_Region' + sr
        init = getattr(config, parname)
        limits = [prior[0], init[0], init[1], prior[1]]
        fixed = prior[2]
        priorshape = prior[3]
        DeltaDec = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # AxialRatio
        AxialRatio = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])
        parname = 'Prior_AxialRatio_Lens' + sl + '_Region' + sr
        prior = getattr(config, parname)
        parname = 'Init_AxialRatio_Lens' + sl + '_Region' + sr
        init = getattr(config, parname)
        limits = [prior[0], init[0], init[1], prior[1]]
        fixed = prior[2]
        priorshape = prior[3]
        AxialRatio = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # PositionAngle
        PositionAngle = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])
        parname = 'Prior_PositionAngle_Lens' + sl + '_Region' + sr
        prior = getattr(config, parname)
        parname = 'Init_PositionAngle_Lens' + sl + '_Region' + sr
        init = getattr(config, parname)
        limits = [prior[0], init[0], init[1], prior[1]]
        fixed = prior[2]
        priorshape = prior[3]
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
        parname = 'Prior_IntrinsicFlux_Source' + ss + '_Region' + sr
        prior = getattr(config, parname)
        parname = 'Init_IntrinsicFlux_Source' + ss + '_Region' + sr
        init = getattr(config, parname)
        limits = [prior[0], init[0], init[1], prior[1]]
        fixed = prior[2]
        priorshape = prior[3]
        IntrinsicFlux = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # DeltaRA
        DeltaRA = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])
        parname = 'Prior_DeltaRA_Source' + ss + '_Region' + sr
        prior = getattr(config, parname)
        parname = 'Init_DeltaRA_Source' + ss + '_Region' + sr
        init = getattr(config, parname)
        limits = [prior[0], init[0], init[1], prior[1]]
        fixed = prior[2]
        priorshape = prior[3]
        DeltaRA = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # DeltaDec
        DeltaDec = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])
        parname = 'Prior_DeltaDec_Source' + ss + '_Region' + sr
        prior = getattr(config, parname)
        parname = 'Init_DeltaDec_Source' + ss + '_Region' + sr
        init = getattr(config, parname)
        limits = [prior[0], init[0], init[1], prior[1]]
        fixed = prior[2]
        priorshape = prior[3]
        DeltaDec = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # Size
        parname = 'Prior_Size_Source' + ss + '_Region' + sr
        prior = getattr(config, parname)
        parname = 'Init_Size_Source' + ss + '_Region' + sr
        init = getattr(config, parname)
        limits = [prior[0], init[0], init[1], prior[1]]
        fixed = prior[2]
        priorshape = prior[3]
        Size = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # AxialRatio
        AxialRatio = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])
        parname = 'Prior_AxialRatio_Source' + ss + '_Region' + sr
        prior = getattr(config, parname)
        parname = 'Init_AxialRatio_Source' + ss + '_Region' + sr
        init = getattr(config, parname)
        limits = [prior[0], init[0], init[1], prior[1]]
        fixed = prior[2]
        priorshape = prior[3]
        AxialRatio = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])

        # PositionAngle
        PositionAngle = dict([('Limits', limits), ('FixedTo', fixed), 
            ('PriorShape', priorshape)])
        parname = 'Prior_PositionAngle_Source' + ss + '_Region' + sr
        prior = getattr(config, parname)
        parname = 'Init_PositionAngle_Source' + ss + '_Region' + sr
        init = getattr(config, parname)
        limits = [prior[0], init[0], init[1], prior[1]]
        fixed = prior[2]
        priorshape = prior[3]
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
