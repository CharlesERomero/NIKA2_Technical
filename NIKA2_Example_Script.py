import Coord_Transforms as CT
import astropy.coordinates as apc
import numpy as np
from astropy import units as u

obj_ra = apc.Angle('12h00m00s')
ojb_dec= apc.Angle('+30d00m00s')
ojbsky = apc.SkyCoord(obj_ra, obj_dec, equinox = 'J2000')

tObs = 400*u.min # I'm currently not using this variable.
#nkotf = CT.nkotf_parameters(XSize=5,YSize=4,PA=0,Tilt=0,Step=10.0,Speed=30.0,
#                            CoordSys='azel',fSamp=20.0)
##### Nico's script:
nkotf = CT.nkotf_parameters(XSize=8.0,YSize=5.0,PA=0,Tilt=0,Step=20.0,Speed=40.0,
                            CoordSys='azel',fSamp=20.0)
#nkotf = CT.nkotf_parameters(XSize=12,YSize=5,PA=0,Tilt=0,Step=10.0,Speed=60.0,
#                            CoordSys='azel',fSamp=20.0)
### I can leave the line below as a commented line for any publicly avaible sample code
pwv = CT.pwv_from_225(0.25) # 20-Oct-2017 reading
pwv = 4.0       # Assume marginal weather.
target='Abell_2443'

##############################################################
Coverage=CT.Multiple_Scans(tObs,a2443sky,nkotf=nkotf,pwv=pwv,elMin=45,doPlot=False)
CT.plot_coverage(Coverage,filename=target+"_above45",target=target)
#                   secObj=g1200p1sky)
CT.plot_visibility(Coverage.date_obs,a2443sky,Coverage,mylabel=target,
                   elMin=45,filename = target+"_Visibility_above45")
CT.hist_pas(Coverage,addname=target)
#CT.make_fits(Coverage,target=target+"_A2443_full_night")

##############################################################
Coverage = CT.Single_Scan(g2sky,nkotf=nkotf,pwv=pwv,elMin=30)
CT.plot_coverage(Coverage,filename=target+"_coverage_1scan",target=target,
                 secObj=g1200p1sky)
CT.plot_visibility(mydate,g2sky,Coverage,mylabel=target,
                   filename = target+"_Visibility_Chart_v3")

##############################################################
target='Abell_2443'
Coverage = CT.Single_Scan(a2443sky,nkotf=nkotf,pwv=pwv,elMin=45)
CT.plot_coverage(Coverage,filename=target+"_coverage_1scan",target=target)
CT.plot_visibility(Coverage.date_obs,a2443sky,Coverage,mylabel=target,
                   elMin=45,filename = target+"_Visibility_Chart_v3")

#goodMin, mydate = CT.find_good_times(elMin=40.0,skyobj=a2443sky)
#CT.plot_visibility(mydate,a2443sky,mylabel="Abell_2443",filename = "A2443_Visibility_Chart_v2")
#TimeAr, ScanX, ScanY = CT.nkotf_scan()
#tStart = CT.create_mytime(mydate,np.min(goodMin))
#RaDecs, ScanAltAz = CT.altaz_SCAN_radec(TimeAr,ScanX,ScanY,tStart,obj=a2443sky)
#CT.plot_radecs(RaDecs,span=600)
#Coverage=None
#Coverage=CT.coverage_map(RaDecs,Coverage=Coverage)
#CT.plot_coverage(Coverage,filename="Abell_2443_coverage_1scan")
#CT.make_fits(Coverage,target="Abell 2443")
CT=reload(CT)
