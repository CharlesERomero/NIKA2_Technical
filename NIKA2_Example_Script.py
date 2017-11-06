import NIKA2_Noise_Estimator as NNE
import astropy.coordinates as apc
import numpy as np
from astropy import units as u
from datetime import datetime
import os

#########################################################################
### Define the coordinates and name of your object:
obj_ra = apc.Angle('12h00m00s')
obj_dec= apc.Angle('+30d00m00s')
skyobj = apc.SkyCoord(obj_ra, obj_dec, equinox = 'J2000')
target='MyObjectName'

#########################################################################
### Define how you are observing your object:
tInt  = 20*u.min    # How much integration time do you want?
elMin = 40          # Minimum elevation (degrees)
nkotf1 = NNE.nkotf_parameters(XSize=6.0,YSize=2.1,PA=20,Tilt=0,Step=5.0,
                             Speed=40.0,CoordSys='azel',fSamp=20.0)
nkotf2 = NNE.nkotf_parameters(XSize=6.0,YSize=2.1,PA=-20,Tilt=0,Step=5.0,
                             Speed=40.0,CoordSys='azel',fSamp=20.0)
nkotf3 = NNE.nkotf_parameters(XSize=6.0,YSize=2.1,PA=40,Tilt=0,Step=5.0,
                             Speed=40.0,CoordSys='azel',fSamp=20.0)
nkotf4 = NNE.nkotf_parameters(XSize=6.0,YSize=2.1,PA=-40,Tilt=0,Step=5.0,
                             Speed=40.0,CoordSys='azel',fSamp=20.0)
scanStrat = [nkotf1,nkotf2,nkotf3,nkotf4]
# Make a list of your scan strategy(ies).

#########################################################################
### Define the conditions (as governed by precipitable water vapor):
### You may give pwv directly as the millimeters of pwv, or
### you may calculate from Tau_225 (the taumeter from the IRAM 30m).
pwv       = NNE.pwv_from_t225(0.3)     # Generally good weather
date_obs  = datetime.strptime('02-07-2019 19:24:08', '%d-%m-%Y %H:%M:%S')
precStart = False

#########################################################################
### Housekeeping:
mydir = os.getcwd()    # If you want to specify a directory, do so here.




#########################################################################
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#########################################################################
###                                                                   ###
###     The user should not need to modify anything below here.       ###
###                                                                   ###
#########################################################################
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#########################################################################




### doPlot keyword is almost deprecated...
Coverage=NNE.Observe_Object(skyobj, nkotf=scanStrat, date=date_obs,
                            precStart=precStart,elMin=elMin,
                            tInt=tInt, pwv=pwv,doPlot=False)

elStr = str(int(elMin))
########################################################################
### The following modules ALL require the variable <Coverage>.
###
### This plots 5 different maps:
NNE.plot_coverage(Coverage,filename=target+"_above"+elStr,target=target,
                  mydir=mydir)

### This makes one plot- the visibility of the object:
NNE.plot_visibility(Coverage.date_obs,skyobj,Coverage,elMin=elMin,
                    mylabel=target,filename = target+"_Visibility_above"+elStr,
                    mydir=mydir)

### This makes one plot- a histogram of parallactic angles:
NNE.hist_pas(Coverage,addname=target,mydir=mydir)

### This writes a fits file with the 5 maps plotted in the plot_coverage
NNE.make_fits(Coverage,target=target+"_full_night",mydir=mydir)

### Shelve the coverage for use later
NNE.shelve_coverage(Coverage,filename='Shelved_Coverage_'+target+'.sav',
                    mydir=mydir)

########################################################################
