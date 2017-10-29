import numpy as np
import astropy.coordinates as apc  # http://www.astropy.org/
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.time import Time as APT
from wcsaxes import WCSAxes        # http://wcsaxes.rtfd.org/
import matplotlib.pyplot as plt
import matplotlib.dates
from matplotlib import colors
from matplotlib import colorbar as cb
import pytz, datetime, os

################################################################################

### Relevant values for IRAM 30-meter telescope:
tmlon = apc.Angle('3:23:55.51 degrees')    # Longitude
tmlat = apc.Angle('37:04:06.29 degrees')   # Latitude
tmalt = 2850.0 * u.m                       # Altitude (meters)
tmpre = 740.0  * u.Pa                      # Pressure (Pascals)
### Create a Earth-Location object:
thirtym = apc.EarthLocation.from_geodetic(tmlon,tmlat, tmalt)

### Should deprecate this. Looks about ready for deprecation.
def_ra = apc.Angle('12h00m00.0s')
def_dec= apc.Angle('37d00m00s')   # Will need to incorporate elevation limits
defsky = apc.SkyCoord(def_ra, def_dec, equinox = 'J2000')

### Where to write files (by default)
cwd    = os.getcwd()  # Get current working directory

################################################################################
### Time-manipulation

def get_dt(tz="Europe/Madrid"):
    local = pytz.timezone (tz)
    my_dt = datetime.datetime.now()
    #   test_dt = my_dt + datetime.timedelta(days=2)
    #   local_dt = local.localize(test_dt, is_dst=None)
    local_dt = local.localize(my_dt, is_dst=None)
    utc_dt = local_dt.astimezone (pytz.utc)

    return utc_dt

def astropyTime_from_datetime(dt):

    myAPT = APT(dt.strftime("%Y-%m-%d %H:%M:%S"))

    return myAPT

################################################################################
### Some simple PWV-to 

def t225_from_pwv(pwv):

    tau = pwv*0.058 + 0.004

    return tau

def pwv_from_t225(tau_225):

    pwv = (tau_225 - 0.004)/0.058

    return pwv

################################################################################
### Astrometry Routines

def get_parallactic_angle(ha, dec, lat=tmlat):
    """
    Calculates the parallactic angle. Many astronomy books will provide info,
    or for info easily retrieved on the internet, look here:
    http://www.gb.nrao.edu/~rcreager/GBTMetrology/140ft/l0058/gbtmemo52/memo52.html

    ---------------
    INPUTS:
    ha:         Hour angle of the source (RA - LST)
    dec:        Declination of the source
    lat:        Latitude of the observing site

    """
    
    #pa = np.arctan(np.cos(lat)*np.sin(az), 
    #               np.sin(lat)*np.cos(el) - np.cos(lat)*np.sin(el)*np.cos(az))
    pa = np.arctan(np.sin(ha)/(np.cos(dec)*np.tan(lat)-np.sin(dec)*np.cos(ha)))

    ### If we needed something beyond +/- pi/2:
    #pa = np.arctan2(np.sin(ha),np.cos(dec)*np.tan(lat)-np.sin(dec)*np.cos(ha))

    return pa

def find_transit(skyobj,myAPT,myloc):

    mylon   = myloc.to_geodetic()[0]
    mylst   = myAPT.sidereal_time('apparent',longitude=mylon)
    ha      = (skyobj.ra - mylst).to("hourangle").value * u.hour
    Transit = myAPT + ha

    return Transit

def find_times_above_el(skyobj,myAPT,myloc,elMin=40):

    Transit  = find_transit(skyobj,myAPT,myloc)
    npts     = 1000
    dTransit = np.linspace(-12, 12, npts)*u.hour
    mytimes  = Transit + dTransit
    ### Can delete the line below, and move to the plotting routine.
    dt_arr   = mytimes.datetime    # Datetime array - for plotting...

    myframe  = apc.AltAz(obstime=mytimes, location=myloc)
    objaltazs = skyobj.transform_to(myframe)
    
    GoodEl = (objaltazs.alt.value > elMin)
    elStart= np.min(mytimes[GoodEl])
    elStop = np.max(mytimes[GoodEl])

    return elStart, elStop, mytimes

########################################################################################
### Scanning-related modules:

class nkotf_parameters:

    def __init__(self,XSize=8,YSize=5,PA=0,Tilt=0,Step=20.0,Speed=40.0,
                 CoordSys="azel",fSamp=20.0):

        self.XSize = XSize
        self.YSize = YSize
        self.PA    = PA
        self.Tilt  = Tilt
        self.Step  = Step
        self.Speed = Speed
        self.CoordSys=CoordSys
        self.fSamp = fSamp

        print "Using the following parameters:"
        print "XSize = ",XSize," arcminutes; YSize = ",YSize," arcminutes"
        print "PA =", PA," degrees; Tilt = ",Tilt," degrees"
        print "Step = ",Step," arcseconds; Speed = ",Speed," arcseconds/second"
        print "in the "+CoordSys+" coordinate system."


def nkotf_scan(XSize=8,YSize=5,PA=0,Tilt=0,Step=20.0,Speed=40.0,
               CoordSys="azel",fSamp=20.0):
    """
    This is meant to create an "X", "Y", and "Time" array for the standard
    Pako script: @nkotf. Currently position angle (PA) and Tilt are *NOT*
    implemented. 

    ----
    INPUTS:

    XSize       - The size of the map along the X direction, in arcminutes
    YSize       - The size of the map along the Y direction, in arcminutes
    PA          - Position Angle
    Tilt        - Tilt of the scans
    Step        - Step in arcseconds
    Speed       - Speed along the X direction (arcseconds per second)
    CoordSys    - The coordinate system used for the scan (either "azel" or "radec")
    fSamp       - A "nominal" sampling frequency (Hz), to indicate how often to mark the
                  trajectory of the telescope. 20 Hertz is the default.
    """
    tTurn  = 3.0    # Seconds
    tTune  = 12.0   # Seconds

    mySign = 1.0
    if (Step < 0):
        mySign = -1

        
    Step   = Step*mySign
    YSpA   = 60.0/Step
    nSub   = int(YSize*YSpA)+1
    
    tTotal = (nSub)*(XSize*60.0/Speed) + (YSize*YSpA)*tTurn + tTune
    tSource= (nSub)*(XSize*60.0/Speed)
    tSS    = XSize*60.0/Speed
    nSSS   = int(tSS*fSamp)

    TimeAr = np.array([]); ScanX = np.array([]); ScanY = np.array([])
    
    for sScan in range(nSub):
    
        ssTimeAr = np.arange(nSSS)/fSamp + tTune + (tTurn+tSS)*sScan
        ssScanX  = np.arange(nSSS)*XSize/float(nSSS)*np.cos(Tilt*u.deg).value -\
                   XSize/2.0 
        ssScanY  = np.zeros(nSSS)*YSize/nSSS*np.sin(Tilt*u.deg).value +\
                   sScan*mySign/YSpA - YSize*mySign/2

        if sScan % 2 == 1:
            ssScanX = np.flipud(ssScanX)
            
        TimeAr = np.append(TimeAr,ssTimeAr)
        ScanX  = np.append(ScanX ,ssScanX)
        ScanY  = np.append(ScanY ,ssScanY)

    myScanX = ScanX*np.cos(PA*u.deg).value + ScanY*np.sin(PA*u.deg).value
    myScanY = ScanY*np.cos(PA*u.deg).value - ScanX*np.sin(PA*u.deg).value

    ScanX = myScanX
    ScanY = myScanY
        
    return TimeAr, ScanX, ScanY

########################################################################################
### A major computational chunk is here: transforming from alt-az to ra-dec.

def altaz_SCAN_radec(TimeAr,ScanX,ScanY,tStart,nkotf,location=thirtym,
                     obj=defsky,doPlot = False):

    ScanTime = TimeAr*u.s + tStart

    objaltaz = obj.transform_to(apc.AltAz(obstime=ScanTime,location=location))
    AltObj = objaltaz.alt
    AzObj  = objaltaz.az

    if nkotf.CoordSys.lower() == 'azel'.lower():
        
        ScanAz   = (ScanX/60.0)*u.deg/np.cos(AltObj) + AzObj
        ScanEl   = (ScanY/60.0)*u.deg + AltObj
        ScanAltAz = apc.AltAz(ScanAz,ScanEl,obstime=ScanTime,location=location)
        RaDecs = ScanAltAz.transform_to(apc.ICRS)

    else:

        # This is at all true! But, this is only used for extinction correction.
        ScanAltAz = objaltaz  
        ScanRa    = (ScanX/60.0)*u.deg/np.cos(obj.dec) + obj.ra
        ScanDec   = (ScanY/60.0)*u.deg + obj.dec
        RaDecs    = apc.SkyCoord(ScanRa, ScanDec, equinox = 'J2000')

    return RaDecs, ScanAltAz

########################################################################################
### Some utility modules:

def get_RaDec_range(RaDec,span=600):

    minRA,maxRA,avgRA = np.min(RaDecs.ra),np.max(RaDecs.ra),np.mean(RaDecs.ra)
    minDEC,maxDEC = np.min(RaDecs.dec.to("deg")),np.max(RaDecs.dec.to("deg"))
    avgDEC = np.mean(RaDecs.dec)
    xrange = [(avgRA -span*u.arcsec).value,(avgRA +span*u.arcsec).value]
    yrange = [(avgDEC-span*u.arcsec).value,(avgDEC+span*u.arcsec).value]

    return xrange,yrange

def int_scalar(x):
    return np.int(x)

int_arr = np.vectorize(int_scalar)

def create_wcs(PixS,avgRA,avgDEC,Xcen,Ycen):

    RAdelt = -PixS.to("deg"); DECdelt = PixS.to("deg")
    w = wcs.WCS(naxis=2)
    # Set up an "gnomic" projection
    # Vector properties may be set with Python lists, or Numpy arrays
    w.wcs.crpix = [Xcen,Ycen] #[avgRA.value,avgDEC.value]
    w.wcs.cdelt = np.array([RAdelt.value,DECdelt.value])
    w.wcs.crval = [avgRA.value,avgDEC.value] #[Xcen,Ycen]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    
    return w

class CovMap:
    
    def __init__(self,in_map, ra0, dec0, pixs,myFOVx,myFOVy,nFOVpix,
                 avgRA,avgDEC,w,radmap,nkotf,date_obs,
                 coordsys='RaDec',notes=None):
        
        self.time   = in_map   # This should stay as a time map (in seconds)
        xsz,ysz     = in_map.shape
        ### OOL FTL IJW TSM
        self.noise1mm  = np.zeros((xsz,ysz))  # A final product, likely in mJy/beam
        self.noise2mm  = np.zeros((xsz,ysz))  # A final product, likely in mJy/beam
        self.weight1mm = np.zeros((xsz,ysz))  # This is initially calculated as a RELATIVE weight
        self.weight2mm = np.zeros((xsz,ysz))  # This is initially calculated as a RELATIVE weight
        self.ra0    = ra0      # Reference pixel [0,0]
        self.dec0   = dec0     # Reference pixel [0,0]
        self.pixs   = pixs     # Pixel Size (units in variable)
        self.myFOVx = myFOVx   # 
        self.myFOVy = myFOVy   #
        self.nFOVpix= nFOVpix  # Number of pixels in the FOV
        self.RAcen  = avgRA    # Should be the center of the map (RA)
        self.DECcen = avgDEC   # Should be the center of the map (Dec)
        self.tint   = 0.0*u.s  # (Total) Integration time (seconds)
        self.w      = w        # WCS structure
        self.radmap = radmap   # Radial map (in whatever units PixS was)
        self.nkotf  = nkotf    # A structure of the OTF parameters used

        self.scannum= np.array([])
        self.scanaz = np.array([])
        self.scanel = np.array([])
        self.scanpwv=  np.array([])
        self.scantau1mm = np.array([])
        self.scantau2mm = np.array([])
        self.scanstart  = np.array([])
        self.scanstop   = np.array([])
        self.scanPA = np.array([])

        ### Some housekeeping info:
        self.date_made = get_dt("Europe/Madrid").strftime("%Y-%m-%d")
        self.date_obs  = date_obs
        if not (notes == None):
            self.notes = notes
        else:
            self.notes = "No comments left."

########################################################################################
### The workhorse routine

def coverage_map(tStart,nkotf,obj,TimeAr, ScanX, ScanY,
                 FoV=6.5,PixS=3.0*u.arcsec,span=None,
                 Coverage=None,fSamp=20.0,pwv=4,myloc=thirtym):
    """
    This is a workhorse script which takes a structure comprised of an array of
    RA positions and an array of Dec positions and creates a map of pixels with
    time spent in/on the *map* pixels. It is simplistic in that it assumes
    COMPLETE, UNIFORM coverage within the FOV. We know this is not true, but a
    more sophisticated treament will have to wait.
    ----
    HISTORY:
    Written 17-Oct-2017 by Charles Romero


    
    Parameters
    ----------
    RaDecs   : A structure from astropy
    Fov      : The Field of View (in arcminutes) of the instrument (NIKA2).
    PixS     : The pixel size used for our coverage map (in arcseconds)
    span     : The (~radial) span (X - Span, X+Span) of the map in arcseconds
    Coverage : A structure including the map and simple astrometric parameters
    fSamp    : Sampling frequency (of our scan)
    Tau      : Zenith tau
    Elevation: Telescope elevation angle (!not altitude in meters!)

    """
    RaDecs, ScanAltAz = altaz_SCAN_radec(TimeAr,ScanX,ScanY,tStart,nkotf,obj=obj)

#    int_arr = np.vectorize(int_scalar)
    Elevation=ScanAltAz.alt
    Azimuth  =ScanAltAz.az
    avgAZ = np.mean(Azimuth)
    avgEL = np.mean(Elevation)
    Tau1mm =opacity_by_band(band="1mm",pwv=pwv)
    Tau2mm =opacity_by_band(band="2mm",pwv=pwv)

    delRa  = (RaDecs.ra - np.roll(RaDecs.ra,1)).to('arcsec').value
    delDec = (RaDecs.dec- np.roll(RaDecs.dec,1)).to('arcsec').value
    negdRa = np.median(delRa[(delRa < 0)])
    posdRa = np.median(delRa[(delRa > 0)])
    negdDec= np.median(delDec[(delRa < 0)])
    posdDec= np.median(delDec[(delRa > 0)])
    scPA   = np.arctan2(negdRa,negdDec)*(u.rad).to("deg")
    scPA2  = np.arctan2(posdRa,posdDec)*(u.rad).to("deg")

    if (scPA > 0) and (scPA < 180):
        myPA = 90.0 - scPA  
    else:
        myPA = 90.0 -scPA2 

    ### Checking my parallactic angle calculations...
    #tMid    = tStart + np.median(TimeAr)*u.s
    #mylon   = myloc.to_geodetic()[0]
    #mylst   = tMid.sidereal_time('apparent',longitude=mylon)
    #ha      = (obj.ra - mylst).to("hourangle")
    #ParAng = get_parallactic_angle(ha,obj.dec)
    #print myPA, ParAng.to("deg")

    scanlen = ((nkotf.XSize/2.0)**2 + (nkotf.YSize/2.0)**2)**0.5
    scanext = scanlen + FoV/1.5
    mybuffer= 1.6

    if span == None:
        span = (scanext*mybuffer*u.arcmin).to('arcsec')

    scanstart = tStart
    scanstop  = tStart + np.max(TimeAr)*u.s
    
    if Coverage == None:
        avgRA ,avgDEC  = np.mean(RaDecs.ra),np.mean(RaDecs.dec)
        nXpix  = int((2*span/PixS).value) ; nYpix = int((2*span/PixS).value)
        #print 'Span = ',span
        #print nXpix,nYpix

        XXarr  = (nXpix/2.0 - np.arange(nXpix))*PixS  # Sky-right coordinates
        YYarr  = (np.arange(nYpix)-nYpix/2.0)*PixS
        XXmap  = np.outer(XXarr, np.zeros(nYpix)+1.0)*u.arcsec
        YYmap  = np.outer(np.zeros(nXpix)+1.0, YYarr)*u.arcsec
        RAmap  = XXmap + avgRA; DECmap = YYmap + avgDEC
        
        #        RAmap  = np.outer(XXarr + AvgRA, np.zeros(nYpix)+1.0)
        #        DECmap = np.outer(np.zeros(nYpix)+1.0,YYarr + AvgDEC)

        RRmap    = (XXmap**2 + YYmap**2)**0.5
        myFOVind = (RRmap < (FoV/2.0)*u.arcmin)
        myFOVx   = XXmap[myFOVind]
        myFOVy   = YYmap[myFOVind]
        nFOVpix  = len(myFOVx)
        refRA    = RAmap[0,0];        refDEC   = DECmap[0,0]
        Xcen = nXpix/2.0; Ycen = nYpix/2.0

        w = create_wcs(PixS,avgRA,avgDEC,Xcen,Ycen)
        Coverage = CovMap(RAmap.value*0.0,refRA,refDEC,PixS,myFOVx,myFOVy,nFOVpix,
                          avgRA,avgDEC,w,RRmap,nkotf,tStart)

        ################################################################
        ### This requires too much memory.
#        RAs2grid = np.outer(RaDecs.ra, np.zeros(nFOVpix)+1.0) # in deg
#        DECs2grid= np.outer(RaDecs.dec, np.zeros(nFOVpix)+1.0)# in deg
#        xFOV2add = np.outer(np.zeros(len(RaDecs.ra)) +1.0, myFOVx)
#        yFOV2add = np.outer(np.zeros(len(RaDecs.dec))+1.0, myFOVy)
#        RAs2grid = (RAs2grid*u.deg - avgRA).to("arcsec")  + xFOV2add 
#        DECs2grid= (DECs2grid*u.deg - avgDEC).to("arcsec")+ yFOV2add

    ExtCorr1mm = np.exp(Tau1mm/np.cos(Elevation)) # Sec(Elevation) approximation for airmass
    ExtCorr2mm = np.exp(Tau2mm/np.cos(Elevation)) # Sec(Elevation) approximation for airmass
    # If we want to do better, we will want to calculate airmass via AltAz somehow??
    # I'm not sure what exists there, but astropy suggests that it can account for atmospheric
    # refraction. Another time...
    if not hasattr(ExtCorr1mm.value, "__len__"):
        ExtCorr1mm = np.zeros(len(RaDecs.ra)) + ExtCorr1mm
        ExtCorr2mm = np.zeros(len(RaDecs.ra)) + ExtCorr2mm

    Coverage.scannum =np.append(Coverage.scannum, len(Coverage.scannum)+1)
    Coverage.scanaz  =np.append(Coverage.scanaz, avgAZ)
    Coverage.scanel  =np.append(Coverage.scanel, avgEL)
    Coverage.scanpwv=  np.append(Coverage.scanpwv,pwv)
    Coverage.scantau1mm = np.append(Coverage.scantau1mm,Tau1mm)
    Coverage.scantau2mm = np.append(Coverage.scantau2mm,Tau2mm)
    Coverage.scanstart = np.append(Coverage.scanstart,scanstart)
    Coverage.scanstop = np.append(Coverage.scanstop,scanstop)
    Coverage.scanPA = np.append(Coverage.scanPA,myPA)

    for i in range(len(RaDecs.ra)):
        RAhitmap  = (Coverage.ra0 - RaDecs.ra[i]).to("arcsec") + Coverage.myFOVx
        DEChitmap = (RaDecs.dec[i]-Coverage.dec0).to("arcsec") + Coverage.myFOVy

        hitmap = [int_arr((RAhitmap/Coverage.pixs).value),
                  int_arr((DEChitmap/Coverage.pixs).value)]
        Coverage.time[hitmap] += 1.0/fSamp
        Coverage.weight1mm[hitmap] += 1.0/(fSamp*(ExtCorr1mm.value[i])**2)
        Coverage.weight2mm[hitmap] += 1.0/(fSamp*(ExtCorr2mm.value[i])**2)
                    
        Coverage.tint += (1.0/fSamp)*u.s

    return Coverage

class nkotf_parameters:

    def __init__(self,XSize=8,YSize=5,PA=0,Tilt=0,Step=20.0,Speed=40.0,
                 CoordSys="azel",fSamp=20.0):

        self.XSize = XSize
        self.YSize = YSize
        self.PA    = PA
        self.Tilt  = Tilt
        self.Step  = Step
        self.Speed = Speed
        self.CoordSys=CoordSys
        self.fSamp = fSamp

        print "Using the following parameters:"
        print "XSize = ",XSize," arcminutes; YSize = ",YSize," arcminutes"
        print "PA =", PA," degrees; Tilt = ",Tilt," degrees"
        print "Step = ",Step," arcseconds; Speed = ",Speed," arcseconds/second"
        print "in the "+CoordSys+" coordinate system."

#############################################################################################
#############################################################################################

def Observe_Object(skyobj, nkotf=None, date=None, precStart=False, elMin=40.0,
                   tInt=24.0*60.0*u.min, pwv=5,doPlot=False, **kwargs):
    
    """
    This is the main program that is called by the astronomer.
    


    ------------------
    INPUTS:

    skyobj:           A sky object, largely just its RA and DEC
    nkotf:            A dictionary/object containing parameters used in an nkotf scan
                      **NOTE** This can also be an array or list of such objects
    date:             A datetime object, in UTC. If none, today is assumed
    precStart:        Set to True if you want your observations to start EXACTLY at the
                      datetime you specify (again, UTC), or if unspecified, then moment that
                      you call this module. Otherwise, scans will start when
                      the object rises above elMin.
    elMin:            Minimum elevation (in degrees, but in Python, just given as a value)
    tInt:             Total *integration* time, in minutes (Does NOT include overhead)!
    pwv:              Precipitable Water Vapor (in mm).
    doPlot:           Depricated???

    ------------------
    NOTES:

    If you wish to observe with many different types of scans, which may depend on elevation,
    opacity, or other variables, this is left to the astronomer to identify when these conditions
    are met, and the astronomer can use the variables <precStart> and <tInt> to constrain the
    times of these observations.

    Currently, there is not a check that one scan is larger than tInt. It will complete, and then
    break from the loop. That is, beware that completion does NOT mean exact compliance with the
    input tInt.

    """
    
    if date == None:
        date = get_dt("Europe/Madrid")      # Currently only made for use with IRAM 30m

    myAPT   = astropyTime_from_datetime(date)    
    elStart, elStop, mytimes = find_times_above_el(skyobj,myAPT,thirtym,elMin=elMin)
    time_per_day = (elStop - elStart).to("minute")

    if precStart == True:
        if elStart < myAPT:
            elStart = myAPT
        else:
            raise Exception("Your precise start date is below the minimum elevation."+
                            "Please change elMin to a lower elevation.")
       
    if nkotf == None:
        nkotf = nkotf_parameters(**kwargs)

    ObsValid = True     #
    Coverage = None     #
    timeInte = 0.0      # Set my integration time to zero
    tStart   = elStart  # The start of a given scan.
    Nloop    = 0
    Nscans   = 0

    ### If anyone is bold enough to predict how pwv changes with time, you could move this
    ### inside the loops.
    
    NEFD1mm, eta1mm, Tau1mm = values_by_band(band="1mm",pwv=pwv)
    NEFD2mm, eta2mm, Tau2mm = values_by_band(band="2mm",pwv=pwv)

    while ObsValid:
        
        for mynkotf in nkotf:
            TimeAr, ScanX, ScanY = nkotf_scan(XSize=mynkotf.XSize,YSize=mynkotf.YSize,PA=mynkotf.PA,
                                              Tilt=mynkotf.Tilt, Step=mynkotf.Step, Speed=mynkotf.Speed,
                                              fSamp=mynkotf.fSamp)
            WallTime   = (np.max(TimeAr) - np.min(TimeAr))*u.s
            Coverage = coverage_map(tStart,mynkotf,skyobj,TimeAr, ScanX, ScanY,
                                    Coverage=Coverage,pwv=pwv)
            tStart  += WallTime
            Nscans  += 1
            ObsValid = (tStart < elStop) & (Coverage.tint < tInt)

            
            if ObsValid:
                tRemainingEl = (elStop - tStart).value * (u.day.to("min"))
                tRemainingInt = (tInt - Coverage.tint).to("min").value
                tRemaining = np.min([tRemainingEl,tRemainingInt])
                StrRem ="{:.2f}".format(tRemaining) + ' minutes'
                StrInt = "{:.2f}".format(Coverage.tint.to("min").value)+" minutes"
                print 'Completed ',Nscans,' scan so far (tInt = '+StrInt+'); '+StrRem+' remaining'
            else:
                StrInt = "{:.2f}".format(Coverage.tint.to("min").value)+" minutes"
                print 'Completed ',Nscans,' scan so far (tInt = '+StrInt+'); this is the last scan.'
                break 

            if doPlot == True:
                plot_radecs(RaDecs)
                import pdb; pdb.set_trace()

        Nloop+=1

    Coverage.weight1mm *= eta1mm/(NEFD1mm**2)
    Coverage.weight2mm *= eta2mm/(NEFD2mm**2)
    nzwt1mm = (Coverage.weight1mm > 0); nzwt2mm = (Coverage.weight2mm > 0)
    Coverage.noise1mm[nzwt1mm] = Coverage.weight1mm[nzwt1mm]**(-0.5)
    Coverage.noise2mm[nzwt2mm] = Coverage.weight2mm[nzwt2mm]**(-0.5)

    return Coverage

#####################################################################################
### NIKA2-specific values:
    
def values_by_band(band="1mm",pwv=4):

    NEFD_1mm = 35.0      # mJy s**1/2
    NEFD_2mm = 10.0      # mJy s**1/2
    eta_1mm  = 0.6       # Fraction of detectors which are used
    eta_2mm  = 0.6       # Fraction of detectors which are used
  
    if band == "1mm":
        NEFD = NEFD_1mm
        eta  = eta_1mm
    if band == "2mm":
        NEFD = NEFD_2mm
        eta  = eta_2mm

    tau = opacity_by_band(band=band,pwv=pwv)

    return NEFD, eta, tau

def opacity_by_band(band="1mm",pwv=4):

    bv_1                  = 0.075   # band1, do not change
    cv_1                  = 0.001   # band1, do not change
    bv_2                  = 0.025   # band2, do not change
    cv_2                  = 0.001   # band2, do not change

    if band == "1mm":
         tau = bv_1*pwv + cv_1
    if band == "2mm":
        tau = bv_2*pwv + cv_2

    return tau

#################################################################################
### Some calculations

def find_radial_noise(Coverage,atR=None,inR=None,profile=False,band="1mm"):

    if band == "1mm":
        nzwt = (Coverage.weight1mm > 0)
        noise = Coverage.noise1mm[nzwt]
        weight= Coverage.weight1mm[nzwt]
    else:
        nzwt = (Coverage.weight2mm > 0)
        noise = Coverage.noise2mm[nzwt]
        weight= Coverage.weight2mm[nzwt]

    myrads  = Coverage.radmap[nzwt]

    index_array = np.argsort(myrads)
    radsorted = myrads[index_array]
    noisesort = noise[index_array]
    weightsort= weight[index_array]
    Rbuffer = Coverage.pixs
    
    if not (atR == None):

#        inRind  = (radsorted < atR)
#        myind   = np.where(radsorted == np.max(radsorted[inRind]))
        inRind  = ((radsorted > atR-Rbuffer) & (radsorted < atR))
        myind   = np.where(radsorted == np.max(radsorted[inRind]))
        if hasattr(myind, "__len__"):
            myind=myind[0][0]
        #mynoise = noisesort[myind];        myweight= weightsort[myind]
        mynoise = np.mean(noisesort[myind])

        return mynoise

    if not (inR == None):

        inRind  = (radsorted < inR)
        myweight= np.mean(weightsort[inRind])
        mynoise = myweight**(-0.5)

        return mynoise

    if (profile == True):

        return radsorted, noisesort, weightsort
        
    
#################################################################################
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#################################################################################
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#################################################################################
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#################################################################################
#################################################################################
##########################                             ##########################
#############                                                       #############
#############                    OUTPUT ROUTINES:                   #############
#############               MAPPING AND PLOTTING ROUTINES           #############
#############                                                       #############
##########################                             ##########################
#################################################################################
#################################################################################
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#################################################################################
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
#################################################################################
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -#
#--  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- #
#################################################################################

def make_fits(Coverage,target="Object",mydir=cwd):

    header = Coverage.w.to_header()

    hdu1 = fits.PrimaryHDU(Coverage.time,header=header)
    hdu1.header.append(("Title","Time Map"))
    hdu1.header.append(("Target",target))
    
    ### Weight Maps:
    hdu2 = fits.ImageHDU(Coverage.weight1mm)
    hdu2.header = header
    hdu2.name = 'Weight_Map_1mm'
    hdu2.header.append(("Title","Weight Map 1mm"));    hdu2.header.append(("Target",target))
    hdu2.header.append(("XTENSION","Second"));         hdu2.header.append(("SIMPLE","T")) 
    hdu2.verify('fix')
    hdu3 = fits.ImageHDU(Coverage.weight2mm)
    hdu3.header = header
    hdu3.name = 'Weight_Map_2mm'
    hdu3.header.append(("Title","Weight Map 1mm"));    hdu3.header.append(("Target",target))
    hdu3.header.append(("XTENSION","Second"));         hdu3.header.append(("SIMPLE","T")) 
    hdu3.verify('fix')
    
    ### Noise Maps:
    hdu4 = fits.ImageHDU(Coverage.noise1mm)
    hdu4.header = header
    hdu4.name = 'Noise_Map_1mm'
    hdu4.header.append(("Title","Noise Map 1mm"));    hdu4.header.append(("Target",target))
    hdu4.header.append(("XTENSION","Second"));        hdu4.header.append(("SIMPLE","T")) 
    hdu4.verify('fix')
    hdu5 = fits.ImageHDU(Coverage.noise2mm)
    hdu5.header = header
    hdu5.name = 'Noise_Map_2mm'
    hdu5.header.append(("Title","Noise Map 1mm"));    hdu5.header.append(("Target",target))
    hdu5.header.append(("XTENSION","Second"));        hdu5.header.append(("SIMPLE","T")) 
    hdu5.verify('fix')

    hdu1.header.add_history("Coverage maps made on "+Coverage.date_made+ ".")
    hdu1.header.add_history("Coverage maps are for observations on "+
                            Coverage.date_obs.datetime.strftime("%Y-%m-%d")+ ".")
    hdulist = fits.HDUList([hdu1,hdu2,hdu3,hdu4,hdu5])
    hdulist.info()
    filename="Coverage_Maps_"+target+".fits"
    fullpath = os.path.join(mydir,filename)
    hdulist.writeto(fullpath,clobber=True,output_verify="exception")

def shelve_coverage(Coverage,filename='Shelved_Coverage.sav',mydir=cwd):

    to_shelve=['Coverage']
    my_shelf = shelve.open(filename,'n') # 'n' for new
    for key in to_shelve:
        if key in cannot_shelve:
            print 'Skipping trying to shelve ',key
        else:
            try:
                my_shelf[key] = globals()[key]
            except TypeError:
                # __builtins__, my_shelf, and imported modules can not be shelved.
                print('ERROR shelving: {0}'.format(key))
    my_shelf.close()


def add_noise_text(Coverage,band="1mm",myfontsize=15,small=False,large=False):
   
    nXpix,nYpix = Coverage.time.shape
    Rdef = 2.0
    if small == True: Rdef=1.0
    if large == True: Rdef=3.0
    noise1mm2amin = find_radial_noise(Coverage,atR=Rdef*u.arcmin,band=band)
    noise1mm4amin = find_radial_noise(Coverage,atR=Rdef*2*u.arcmin,band=band)
    noise1mm6amin = find_radial_noise(Coverage,atR=Rdef*3*u.arcmin,band=band)
    Rone,Rtwo,Rthr = Rdef,Rdef*2,Rdef*3
    Rons,Rtws,Rths = int(Rone),int(Rtwo),int(Rthr)
    StrOns,StrTws,StrThs = str(Rons),str(Rtws),str(Rths)
    plt.text(10.0*nXpix/400,10.0*nYpix/400,"Noise in mJy/beam",color='blue',fontsize=myfontsize)
    plt.text(10.0*nXpix/400,30.0*nYpix/400,r'$\sigma_{'+StrOns+'^{\prime}} = $ '+"{:.3f}".format(noise1mm2amin),
             color='red',fontsize=myfontsize)
    plt.text(10.0*nXpix/400,50.0*nYpix/400,r'$\sigma_{'+StrTws+'^{\prime}} = $ '+"{:.3f}".format(noise1mm4amin),
             color='orange',fontsize=myfontsize)
    plt.text(10.0*nXpix/400,70.0*nYpix/400,r'$\sigma_{'+StrThs+'^{\prime}} = $ '+"{:.3f}".format(noise1mm6amin),
             color='green',fontsize=myfontsize)
    plt.contour(Coverage.radmap.to("arcmin").value, [Rone,Rtwo,Rthr],
                colors=('red','orange','green'))

    inthours = "{:.2f}".format(Coverage.tint.to("hour").value)
    plt.text(10.0*nXpix/400,370.0*nYpix/400,r'$t_{int} = $'+inthours+' hours',fontsize=myfontsize)

def add_scan_text(Coverage,band="1mm",myfontsize=15):

    nXpix,nYpix = Coverage.time.shape
    AllAvgEl  = np.mean(Coverage.scanel.value)
    AllAvgPWV = np.mean(Coverage.scanpwv)
    AllAvgTau1= np.mean(Coverage.scantau1mm)
    AllAvgTau2= np.mean(Coverage.scantau2mm)
    AllAvgEC1 = np.mean(np.exp(Coverage.scantau1mm/np.cos(Coverage.scanel.value*(u.deg.to('rad')))))
    AllAvgEC2 = np.mean(np.exp(Coverage.scantau2mm/np.cos(Coverage.scanel.value*(u.deg.to('rad')))))
    
    plt.text(220.0*nXpix/400,10.0*nYpix/400,r'$\langle $Elev$ \rangle = $'+
             "{:.2f}".format(AllAvgEl)+' deg.',fontsize=myfontsize)

    if band == "1mm":
        plt.text(260.0*nXpix/400,350.0*nYpix/400,
                 r'$\langle \tau_{1mm} \rangle = $'+"{:.2f}".format(AllAvgTau1),fontsize=myfontsize)
        plt.text(10.0*nXpix/400,350.0*nYpix/400,
                 r'$\langle $Ext. Corr$_{1mm} \rangle= $'+"{:.2f}".format(AllAvgEC1),
                 fontsize=myfontsize)
    else:
        plt.text(260.0*nXpix/400,350.0*nYpix/400,
                 r'$\langle \tau_{2mm} \rangle = $'+"{:.2f}".format(AllAvgTau2),fontsize=myfontsize)
        plt.text(10.0*nXpix/400,350.0*nYpix/400,
                 r'$\langle $Ext. Corr$_{2mm} \rangle= $'+"{:.2f}".format(AllAvgEC2),
                 fontsize=myfontsize)

    plt.text(230.0*nXpix/400,370.0*nYpix/400,
             r'$\langle $PWV$ \rangle = $'+"{:.2f}".format(AllAvgPWV)+' mm',fontsize=myfontsize)


def plot_skyCoord(myax,Coverage,skyobj,mycolor='purple',mylabel='Other'):


#    xcen,ycen = Coverage.w.wcs_world2pix(g2_ra,g2_dec,0)
#    xxce,yyce = Coverage.w.wcs_world2pix(skyobj.ra,skyobj.dec,0)
#    delx = xcen - xxce
#    dely = ycen - yyce
    delra  = ((skyobj.ra.value  -  Coverage.w.wcs.crval[0])*u.deg).to('arcsec')
    deldec = ((skyobj.dec.value - Coverage.w.wcs.crval[1])*u.deg).to('arcsec')
    delxx  = -delra*np.cos(Coverage.w.wcs.crval[1]*u.deg)/Coverage.pixs
    delyy  = deldec/Coverage.pixs

    myxx = Coverage.w.wcs.crpix[0] + delxx.value
    myyy = Coverage.w.wcs.crpix[1] + delyy.value

    #myrr = ((delxx.value)**2 + (delyy.value)**2)**0.5
    
    myax.plot([myxx],[myyy],'x',color='white',ms=5,mew=2)
    myax.plot([myxx],[myyy],'x',color=mycolor,ms=3,mew=1,label=mylabel)

def hist_pas(Coverage,dpi=200,filename="Parallactic_Angle_Histogram_",
             addname="Object",myfontsize=15,format='png',mydir=cwd):

    fig = plt.figure(dpi=dpi,figsize=(8,8))
    ax = fig.add_axes([0.10, 0.3, 0.85, 0.65])
    ax1 = fig.add_axes([0.10, 0.10, 0.85, 0.10])
    # N is the count in each bin, bins is the lower-limit of the bin
    N, bins, patches = ax.hist(Coverage.scanPA,bins='auto')
    binsize = np.median(bins - np.roll(bins,1))
    binEl  = np.array([]);  paStart= np.array([]); paStop = np.array([])

    for pa in bins:
        gi = (Coverage.scanPA >= pa) & (Coverage.scanPA < pa+binsize)
        myEls = Coverage.scanel[gi].value
        if any(gi):
            binEl  = np.append(binEl,np.mean(myEls))
            paStart = np.append(paStart,np.min(Coverage.scanstart[gi]))
            paStop  = np.append(paStop,np.max(Coverage.scanstop[gi]))
        
    cmap=plt.cm.spectral
    mycolors=np.array([])
    norm = colors.Normalize(binEl.min()*0.98,binEl.max()*1.02)
    for myEl, mypatch in zip(binEl,patches):
        color = cmap(norm(myEl))
        mypatch.set_facecolor(color)
        mycolors = np.append(mycolors,color)

    for pa,mystart in zip(bins,paStart):
        nmax = np.max(N)
        yy = (pa - np.min(bins))*nmax/(np.max(bins) - np.min(bins))
        mytime =  mystart.datetime
        myutc  = mytime.strftime('%H:%M')+' UTC'
        ax.text(pa,yy/2.0+4,myutc,color='black',fontsize=myfontsize,rotation=-90)

    startdt = datetime.datetime.strptime(Coverage.scanstart[0].value, '%Y-%m-%d %H:%M:%S.%f')
    startday= startdt.strftime('%Y-%m-%d')
    myxloc = np.min(bins)*0.6 + np.max(bins)*0.4
    
    ax.text(myxloc,0.9*np.max(N), startday,color='black',fontsize=myfontsize)
    
    units="degrees"
    cb1 = cb.ColorbarBase(ax1, cmap=cmap,norm=norm,
                          orientation='horizontal')
    ax.set_title("Histogram of Parallactic Angles on "+addname,fontsize=myfontsize)
    ax.set_xlabel("Parallactic Angle (degrees)",fontsize=myfontsize)
    ax.set_ylabel("Number of Scans",fontsize=myfontsize)
    ax1.set_xlabel("Average Elevation",fontsize=myfontsize)
    fullbase = os.path.join(mydir,filename)
    fulleps = fullbase+addname+'.eps'; fullpng = fullbase+addname+'.png'

    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')
    plt.clf()

#############################################################################
    
def ind_plots_cov(Coverage,map,filename="NIKA2_Coverage_map",target="Object",
                  myfontsize=15,mytitle="my map",addname="_quantity",
                  units='(units)',band="1mm",cblim=False,addtext=False,dpi=200,
                  secObj=None,thiObj=None,fouObj=None,format='png',mydir=cwd):

    if band == "1mm":
        nzwt = (Coverage.weight1mm > 0)
        mymin = np.min(Coverage.noise1mm[nzwt])
        gwmax = np.max(Coverage.noise1mm[nzwt])*0.5
        mymax = mymin*10.0
    else:
        nzwt = (Coverage.weight2mm > 0)
        mymin = np.min(Coverage.noise2mm[nzwt])
        gwmax = np.max(Coverage.noise2mm[nzwt])*0.5
        mymax = mymin*10.0

    small = False
    large = False
    if (Coverage.nkotf.XSize < 8) and (Coverage.nkotf.YSize < 8):
        small = True
        
    fig = plt.figure(dpi=dpi,figsize=(8,8)); axpos=[0.2, 0.2, 0.7, 0.7]
    ax = WCSAxes(fig, axpos, wcs=Coverage.w);fig.add_axes(ax)
#    cax = ax.imshow(map,interpolation='none',
#                    norm=colors.LogNorm(vmin=mymin,vmax=mymax),cmap='bwr')
    if cblim == False:
        cax = ax.imshow(map,interpolation='none',cmap='bwr',origin='lower')
        plt.contour(map, [0],colors=('black'),linewidths=3)        
    else:
        cax = ax.imshow(map,interpolation='none',origin='lower',
                        norm=colors.LogNorm(vmin=mymin,vmax=mymax),cmap='bwr')
    plt.title(mytitle+target,fontsize=myfontsize*1.2,y=1.08)
    strxs = "{:.1f}".format(float(Coverage.nkotf.XSize))
    strys = "{:.1f}".format(float(Coverage.nkotf.YSize))
    strpa = "{:.1f}".format(float(Coverage.nkotf.PA))
    strti ="{:.1f}".format(float(Coverage.nkotf.Tilt))
    strst = "{:.1f}".format(float(Coverage.nkotf.Step))
    strsp="{:.1f}".format(float(Coverage.nkotf.Speed))
    strofparams = "XSize = "+strxs+", "+"YSize = "+strys+", "+\
                  "PA ="+strpa+", "+"Tilt ="+strti+", "+\
                "Step ="+strst+", "+"Speed ="+strsp+", "
    strnkotf = "@nkotf "+strxs+' '+strys+' '+strpa+' '+strti+' '+strst+' '+\
               strsp+' '+Coverage.nkotf.CoordSys
    plt.suptitle(strnkotf,x=0.48,y=0.87,fontsize=myfontsize,color='blue')
    fullbase = os.path.join(mydir,filename)
    fulleps = fullbase+addname+'.eps'; fullpng = fullbase+addname+'.png'
    cbar = fig.colorbar(cax) ; cbar.set_label(units,fontsize=myfontsize)
    ra = ax.coords[0]; dec = ax.coords[1]
    ra.set_major_formatter('hh:mm:ss.s');dec.set_major_formatter('dd:mm:ss.s')
    ra.set_axislabel("RA (J2000)",fontsize=myfontsize)
    dec.set_axislabel("Dec (J2000)",fontsize=myfontsize)
    if addtext == True:
        add_noise_text(Coverage,band=band,small=small,large=large)
        add_scan_text(Coverage,band=band)

    if secObj != None:
        plot_skyCoord(ax,Coverage,secObj,mylabel='G1200.1')
        
    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')
    plt.clf()   
    
def plot_coverage(Coverage,filename="NIKA2_Coverage_map",target="Object",
                  secObj=None,thiObj=None,fouObj=None,format='png',mydir=cwd):

    myfontsize=15
    ind_plots_cov(Coverage,Coverage.time,filename=filename,target=target,
                  mytitle="Time Map; ",addname="_time",units='seconds',
                  myfontsize=myfontsize,secObj=secObj,thiObj=thiObj,
                  fouObj=fouObj,format=format,mydir=mydir)
    
    ind_plots_cov(Coverage,Coverage.weight1mm,filename=filename,target=target,
                  mytitle="Weight Map, 1mm; ",addname="_weight1mm",units="mJy/beam $^{-2}$",
                  myfontsize=myfontsize,band="1mm",secObj=secObj,thiObj=thiObj,
                  fouObj=fouObj,format=format,mydir=mydir)

    ind_plots_cov(Coverage,Coverage.weight2mm,filename=filename,target=target,
                  mytitle="Weight Map, 2mm; ",addname="_weight2mm",units="mJy/beam $^{-2}$",
                  myfontsize=myfontsize,band="2mm",secObj=secObj,thiObj=thiObj,
                  fouObj=fouObj,format=format,mydir=mydir)

    ind_plots_cov(Coverage,Coverage.noise1mm,filename=filename,target=target,
                  mytitle="Noise Map, 1mm; ",addname="_noise1mm",units="mJy/beam",
                  myfontsize=myfontsize,band="1mm",cblim=True,addtext=True,secObj=secObj,thiObj=thiObj,
                  fouObj=fouObj,format=format,mydir=mydir)

    ind_plots_cov(Coverage,Coverage.noise2mm,filename=filename,target=target,
                  mytitle="Noise Map, 2mm; ",addname="_noise2mm",units="mJy/beam",
                  myfontsize=myfontsize,band="2mm",cblim=True,addtext=True,secObj=secObj,thiObj=thiObj,
                  fouObj=fouObj,format=format,mydir=mydir)

def plot_visibility(mydate,skyobj,Coverage=0,elMin=40,mylabel="Target",
                    dpi=200,filename = "Visibility_Chart",format='png',
                    myloc=thirtym,mydir=cwd):

    elStart, elStop, mytimes = find_times_above_el(skyobj,mydate,myloc,elMin=elMin)
    
    myframe  = apc.AltAz(obstime=mytimes, location=myloc)
    sunaltazs = apc.get_sun(mytimes).transform_to(myframe)
    #moonaltazs= apc.get_moon(mytimes).transform_to(myframe)
    objaltazs = skyobj.transform_to(myframe)

    #############################

    bt30 = (np.max(mytimes[(objaltazs.alt.value > 30.0)]) - \
           np.min(mytimes[(objaltazs.alt.value > 30.0)])).to("hour").value
    bt40 = (np.max(mytimes[(objaltazs.alt.value > 40.0)]) - \
           np.min(mytimes[(objaltazs.alt.value > 40.0)])).to("hour").value
    bt50 = (np.max(mytimes[(objaltazs.alt.value > 50.0)]) - \
           np.min(mytimes[(objaltazs.alt.value > 50.0)])).to("hour").value
    bteM = (np.max(mytimes[(objaltazs.alt.value > elMin)]) - \
           np.min(mytimes[(objaltazs.alt.value > elMin)])).to("hour").value

    plt.figure(1,dpi=dpi,figsize=(8,8));    plt.clf();    fig1,ax1 = plt.subplots()

    date_arr = mytimes.datetime       # Convert to datetime array, for plotting
    
    ax1.fill_between(date_arr, 0, 90,
                     sunaltazs.alt < -0*u.deg, color='0.5', zorder=0)
    ax1.fill_between(date_arr, 0, 90,
                     sunaltazs.alt < -18*u.deg,color='k', zorder=0)
    ax1.plot_date(date_arr, sunaltazs.alt.value, '-',color='r',lw=3, label='Sun')
#    ax1.plot_date(date_arr, moonaltazs.alt.value,'-',color='b',lw=3, label='Moon')
    ax1.plot_date(date_arr, objaltazs.alt.value, '-',color='g',lw=5, label=mylabel)
    myxlim = ax1.get_xlim(); ax1.set_xlim(myxlim)
    ax1.plot(myxlim,[30,30],'--',color ='0.75',label="{:.1f}".format(bt30)+" hrs")
    ax1.plot(myxlim,[40,40],'--',color ='0.5' ,label="{:.1f}".format(bt40)+" hrs")
    ax1.plot(myxlim,[50,50],'--',color ='0.25',label="{:.1f}".format(bt50)+" hrs")
    myStart = elStart.datetime; myStop = elStop.datetime
    
    myylim = ax1.get_ylim(); ax1.set_ylim(myylim)
    ax1.plot(myxlim,[elMin,elMin],'--',color ='b',
             label="{:.1f}".format(bteM)+" hrs above "+str(int(elMin)))
    ax1.plot_date([myStart,myStart],[0,90],'--',color='b')
    ax1.plot_date([myStop,myStop],[0,90],'--',color='b')

    
    gi = np.array([])
    if isinstance(Coverage,NNE.CovMap):
        for start,stop in zip(Coverage.scanstart,Coverage.scanstop):
            mgi = np.where(((mytimes > start) & (mytimes < stop)) == True)
            gi = np.append(gi,mgi)

        gi = int_arr(gi)
        ax1.plot_date(date_arr[gi], objaltazs.alt.value[gi],color='orange',ms=2,
                      label="Obs. ("+"{:.1f}".format(Coverage.tint.to('hour').value)+' hrs)')

        
#    plt.colorbar().set_label('Azimuth [deg]')
    plt.legend(loc='upper left')
    plt.title("Visibility of "+mylabel)
    plt.ylim(0, 90)
    plt.xlabel('UTC (MM-DD HH)')
    plt.ylabel('Altitude [deg]')  
    plt.gcf().autofmt_xdate()  # Hopefully make the x-axis (dates) look better
    fullbase = os.path.join(mydir,filename)
    fulleps = fullbase+'.eps'; fullpng = fullbase+'.png'
    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')

def plot_radecs(RaDecs,span=600,format='png',mydir=cwd):
    """
    This module plots the right asciension and declination of a (raster) scan
    ----
    INPUTS:
    RaDecs : is structure array. RaDecs.ra is an array of quantities (should be
    in degrees) of the right ascensions; RaDecs.dec is the same for declination.

    span   : is given in arcseconds and determines the plot range from the center as
    mid-span to mid+span
    """
    
    plt.figure(figsize=(10,10))
    plt.plot(RaDecs.ra.value,RaDecs.dec.value,'.')
#    cwd    = os.getcwd()  # Get current working directory
#    minRA,maxRA = np.min(RaDecs.ra.value),np.max(RaDecs.ra.value)
#    minDEC,maxDEC = np.min(RaDecs.dec.value),np.max(RaDecs.dec.value)
    minRA,maxRA,avgRA = np.min(RaDecs.ra),np.max(RaDecs.ra),np.mean(RaDecs.ra)
    minDEC,maxDEC = np.min(RaDecs.dec.to("deg")),np.max(RaDecs.dec.to("deg"))
    avgDEC = np.mean(RaDecs.dec)
    xrange = [(avgRA -span*u.arcsec).value,(avgRA +span*u.arcsec).value]
    yrange = [(avgDEC-span*u.arcsec).value,(avgDEC+span*u.arcsec).value]
    plt.xlim(xrange);    plt.ylim(yrange)
    plt.title("Example; One Alt-Az scan + tracking to Ra-Dec")
    
    actSpan = np.max([(maxRA - minRA).value, (maxDEC - minDEC).value])
    if 2*span < actSpan*3600 + 6.25*60:
        print span, actSpan*3600 + 6.25*60
        import pdb; pdb.set_trace()
        #raise Exception("Your plotting range (span) is too small")
    ### Maybe there is further optimizing to do here?
    
    numberofTicks = 5

    beststep = int(span/(numberofTicks+1))
    lowRA = (int(minRA.to("arcsec").value/beststep)+1)*beststep*u.arcsec
    highRA= int(maxRA.to("arcsec").value/beststep)*beststep*u.arcsec
    mminRA = lowRA.to("deg"); mmaxRA = highRA.to("deg")

#    ticksLocationsRA = linspace(minRA, maxRA, numberOfTicks)
#    ticksLabelsRA = degrees_to_hhmmss(ticksLocatoinsRA)
#    xticks(ticksLocationsRA, ticksLabelsRA)

    lowDEC = (int(minDEC.to("arcsec").value/beststep)+1)*beststep*u.arcsec
    highDEC= int(maxDEC.to("arcsec").value/beststep)*beststep*u.arcsec
    mminDEC = lowDEC.to("deg"); mmaxDEC = highDEC.to("deg")
#    ticksLocationsDEC = linspace(minDEC, maxDEC, numberOfTicks)
#    ticksLabelsDEC = degrees_to_ddmmss(ticksLocatoinsDEC)
#    yticks(ticksLocationsDEC, ticksLabelsDEC)  
#    plt.yticks(ticksLocationsDEC,ticksLabelsDEC)
#    plt.xticks(ticksLocationsRA ,ticksLabelsRA)
    
    filename = "RaDec_map_v2";fullbase = os.path.join(mydir,filename)
    fulleps = fullbase+'.eps'; fullpng = fullbase+'.png'
    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')

def plot_altaz(ScanAz,ScanEl,TestRaDec,format='png',mydir=cwd):

    plt.figure()
    plt.plot(ScanAz.value,ScanEl.value,'.')
    filename = "AltAz_map_v2";fullbase = os.path.join(mydir,filename)
    fulleps = fullbase+'.eps'; fullpng = fullbase+'.png'
    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')
        
    plt.figure()
    plt.plot(TestRaDec.ra.value,TestRaDec.dec.value,'.')
    filename = "TestRaDec_map_v2";fullbase = os.path.join(mydir,filename)
    fulleps = fullbase+'.eps'; fullpng = fullbase+'.png'
    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')



