# NIKA2_Technical

This repository was developed with Python 2.7.12 and IPython 2.4.1
The repository works on Python versions XXX and IPython YYY

Third party packages used are:
(1) Astropy: http://www.astropy.org/
(2) WCSaxes: http://wcsaxes.rtfd.org/

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

It is necessary to edit variables above L50 (Release Oct 29, 2017) in the
file <NIKA2_Example_Script.py> (unless, by chance, you want some of the
same parameters).

