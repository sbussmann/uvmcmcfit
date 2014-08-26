Input Data Format
=================

``uvmcmcfit`` uses visibilities obtained from an interferometer as the primary
input data.  The format must be either *uvfits* (preferred) or CASA *ms*.  

Which format should I use?
--------------------------

Use the CASA *ms* format if you are only willing to run ``uvmcmcfit`` inside
the CASA shell.  Doing this will prevent you from taking advantage of the
parallel processing features included in ``uvmcmcfit`` and therefore is likely
to slow down 
