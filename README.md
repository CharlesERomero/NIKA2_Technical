# NIKA2_Technical

This repository was developed with Python 2.7.12 and IPython 2.4.1
The repository has been tested on a variety of Python 2.7 and IPython
2.4 versions. Functionality beyond these versions is not guaranteed. 

Third party packages used are:
(1) Astropy: http://www.astropy.org/
(2) FITS_tools: http://fits-tools.readthedocs.io/en/latest/index.html

Both should be installable with Anaconda or PIP.

Upon using the routine, Astropy will download information from
the US Naval Observatory (for AltAz / RaDec transformations).

####################################################################

To use the package, you can simply run

   > python NIKA2_Example_Script

However, you may find it more useful to use ipython, and copy and paste
lines from NIKA2_Example_Script:

   > ipython

   In [1]: import NIKA2_Noise_Estimator as NNE

   In [2]: import astropy.coordinates as apc

   In [3]: import numpy as np

   In [4]: from astropy import units as u

   In [5]: from datetime import datetime

   In [6]: 

Etc...

It is necessary to edit variables above L45 (Release Jan 20, 2018) in the
file <NIKA2_Example_Script.py> (unless, by chance, you want some of the
same parameters).

