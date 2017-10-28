import numpy as np
import astropy.coordinates as apc
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.time import Time as APT
import time
from wcsaxes import WCSAxes  # http://wcsaxes.rtfd.org/
import matplotlib.pyplot as plt
import matplotlib.dates
from matplotlib import colors
from matplotlib import colorbar as cb
import os
import datetime
#from astropy.wcs.utils import celestial_frame_to_wcs

tmlon = apc.Angle('3:23:55.51 degrees')
tmlat = apc.Angle('37:04:06.29 degrees')
tmalt = 2850.0 * u.m    
thirtym = apc.EarthLocation.from_geodetic(tmlon,tmlat, tmalt)
utcoffset = 2.0*u.hour
mydate  = time.strftime("%Y-%m-%d")
mytoday = mydate

### Should deprecate this. Looks about ready for deprecation.
def_ra = apc.Angle('12h00m00.0s')
def_dec= apc.Angle('37d00m00s')   # Will need to incorporate elevation limits
defsky = apc.SkyCoord(def_ra, def_dec, equinox = 'J2000')

ra=apc.Angle('0h0m0s')
dec= apc.Angle('0d0m0s')
equinox='J2000'
cwd    = os.getcwd()  # Get current working directory

def 225_from_pwv(pwv):

    tau = pwv*0.058 + 0.004

    return tau

def pwv_from_225(tau_225):

    pwv = (tau_225 - 0.004)/0.058

    return pwv

def create_mytime(mydate,myminute,utcoffset=2.0*u.hour):
    
    hours = int(myminute/60.0)
    addday = int(hours/24)
    hours  = hours % 24
    minutes = int(myminute) - hours*60 - addday*24*60
    seconds = int((minutes*60)%60)
    hms = "%02d:%02d:%02d" % (hours,minutes,seconds)
    mytime = APT(mydate+' '+hms) - utcoffset
    mytime = mytime + addday*u.day
    
    return mytime

def find_transit(skyobj,date,utcoffset=2.0*u.hour):

    ntimes=48 # Every half hour
    midnight = create_mytime(date,0,utcoffset)
    delta_midnight = np.linspace(0, 24, ntimes)*u.hour
    mytimes = midnight + delta_midnight
    objaltaz = skyobj.transform_to(apc.AltAz(obstime=mytimes,location=thirtym))
    myels = (objaltaz.alt).to('deg').value

    bestind = np.where(objaltaz.alt == np.max(objaltaz.alt))
    besttime= mytimes[bestind]

    return besttime

def find_good_times(elMin=40,skyobj=None,date=None):

    ra = def_ra   # Will want to delete this...
    dec=def_dec   # Will want to delete this...

    if skyobj == None:
        mysky = apc.SkyCoord(ra, dec, equinox = equinox)
    else:
        mysky = skyobj
        
    if (date == None):
        mydate  = time.strftime("%Y-%m-%d")
    else:
        mydate = date

    scandur= 5.0 # minutes
    ntimes = int(24*60.0/scandur)
    myels = np.zeros(ntimes)

    bestdate = find_transit(skyobj,mydate)
    
#    midnight = create_mytime(mydate,0,utcoffset=utcoffset)
    delta_besttime = np.linspace(-12, 12, ntimes)*u.hour
    mytimes = bestdate + delta_besttime
    objaltaz = mysky.transform_to(apc.AltAz(obstime=mytimes,location=thirtym))
    myels = (objaltaz.alt).to('deg').value
    
######################################################################
#    for sCount in range(ntimes):
#        myminute = sCount*scandur
#        mytime = create_mytime(mydate,myminute,utcoffset=utcoffset)
#        
#        defaltaz = mysky.transform_to(apc.AltAz(obstime=mytime,location=thirtym))
#
#        myels[sCount] = ((defaltaz.alt).to('deg')).value
######################################################################

    goodEl = (myels > elMin)
    myGoodEl = myels[goodEl]
    myGoodT = len(myGoodEl)*scandur

    goodMin = np.arange(ntimes)[goodEl] * scandur
    # Maybe we want to track the object from its rise to its set.
    goodMin = resolve_midnight(goodMin)

#    import pdb;pdb.set_trace()
    goodStart = mytimes[np.min(np.where(goodEl == True))]

    ### I think I want to return myGoodT instead of goodMin
    return goodMin, mydate, goodStart

def resolve_midnight(goodMin):

    myshift = goodMin - np.roll(goodMin,1)
    mystep  = np.median(myshift)
    mybreaks= (myshift != mystep)
    mybrmin = goodMin[mybreaks]
    goodMorn= (goodMin < np.max(mybrmin))
    mynewMin= np.append(goodMin[~goodMorn],goodMin[goodMorn]+24*60.0)

    return mynewMin

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
            ssScanX = np.flip(ssScanX,0)
            
        TimeAr = np.append(TimeAr,ssTimeAr)
        ScanX  = np.append(ScanX ,ssScanX)
        ScanY  = np.append(ScanY ,ssScanY)

    myScanX = ScanX*np.cos(PA*u.deg).value + ScanY*np.sin(PA*u.deg).value
    myScanY = ScanY*np.cos(PA*u.deg).value - ScanX*np.sin(PA*u.deg).value

    ScanX = myScanX
    ScanY = myScanY
        
    return TimeAr, ScanX, ScanY

def altaz_SCAN_radec(TimeAr,ScanX,ScanY,tStart,location=thirtym,
                     obj=defsky,doPlot = False):

    ScanTime = TimeAr*u.s + tStart

    objaltaz = obj.transform_to(apc.AltAz(obstime=ScanTime,location=location))
    AltObj = objaltaz.alt
    AzObj  = objaltaz.az
    # Cross-Elevation to actual Azimuth offsets + Azimuth_0
    ScanAz   = (ScanX/60.0)*u.deg/np.cos(AltObj) + AzObj
    # And actual elevations of the scan
    ScanEl   = (ScanY/60.0)*u.deg + AltObj

    ScanAltAz = apc.AltAz(ScanAz,ScanEl,obstime=ScanTime,location=location)
    
    #######################################################################
    ### TESTING

    TestAz    = AzObj + ScanX*0*u.deg
    TestEl    = AltObj+ ScanY*0*u.deg
    TestAltAz = apc.AltAz(TestAz,TestEl,obstime=ScanTime,location=location)
    TestRaDec = TestAltAz.transform_to(apc.ICRS)

#    print TestRaDec.ra.value[:250]
#    print TestRaDec.dec.value[:250]
#    import pdb; pdb.set_trace()

    if doPlot == True:
        plt.figure()
        plt.plot(ScanAz.value,ScanEl.value,'.')
        filename = "AltAz_map_v2";fullbase = os.path.join(cwd,filename)
        fulleps = fullbase+'.eps'; fullpng = fullbase+'.png'
        plt.savefig(fulleps,format='eps'); plt.savefig(fullpng,format='png')
        
        plt.figure()
        plt.plot(TestRaDec.ra.value,TestRaDec.dec.value,'.')
        filename = "TestRaDec_map_v2";fullbase = os.path.join(cwd,filename)
        fulleps = fullbase+'.eps'; fullpng = fullbase+'.png'
        plt.savefig(fulleps,format='eps'); plt.savefig(fullpng,format='png')

    RaDecs = ScanAltAz.transform_to(apc.ICRS)

    return RaDecs, ScanAltAz

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
#    w = wcs.wcs.utils.celestial_frame_to_wcs(wcs.ICRS(), projection='TAN')
    w = wcs.WCS(naxis=2)
    # Set up an "gnomic" projection
    # Vector properties may be set with Python lists, or Numpy arrays
#    import pdb; pdb.set_trace()
    w.wcs.crpix = [Xcen,Ycen] #[avgRA.value,avgDEC.value]
    w.wcs.cdelt = np.array([RAdelt.value,DECdelt.value])
    w.wcs.crval = [avgRA.value,avgDEC.value] #[Xcen,Ycen]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
#    w.wcs.set_pv([(2,1,45.0)]) # phip, phi0, and theta0 ???
    
    return w

class CovMap:
    
    def __init__(self,in_map, ra0, dec0, pixs,myFOVx,myFOVy,nFOVpix,
                 avgRA,avgDEC,w,radmap,nkotf,
                 coordsys='RaDec',date_made=None,date_obs=None,notes=None):
        
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
        self.scanstart = np.array([])
        self.scanstop  = np.array([])
        self.scanPA = np.array([])

        ### Some housekeeping info:
        if date_made == None:
            date_made = mytoday
        if date_obs == None:
            date_obs = mydate
        self.date_made = date_made
        self.date_obs  = date_obs
        if not (notes == None):
            self.notes = notes
        else:
            self.notes = "No comments left."


def coverage_map(tStart,nkotf,obj,TimeAr, ScanX, ScanY,
                 FoV=6.5,PixS=3.0*u.arcsec,span=None,
                 Coverage=None,fSamp=20.0,pwv=4):
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
    RaDecs, ScanAltAz = altaz_SCAN_radec(TimeAr,ScanX,ScanY,tStart,obj=obj)

#    int_arr = np.vectorize(int_scalar)
    Elevation=ScanAltAz.alt
    Azimuth  =ScanAltAz.az
    avgAZ = np.mean(Azimuth)
    avgEL = np.mean(Elevation)
    Tau1mm =opacity_by_band(band="1mm",pwv=pwv)
    Tau2mm =opacity_by_band(band="2mm",pwv=pwv)

    delRa = (RaDecs.ra - np.roll(RaDecs.ra,1)).to('arcsec').value
    delDec= (RaDecs.dec- np.roll(RaDecs.dec,1)).to('arcsec').value
    negdRa= np.median(delRa[(delRa < 0)])
    posdRa= np.median(delRa[(delRa > 0)])
    negdDec=np.median(delDec[(delRa < 0)])
    posdDec=np.median(delDec[(delRa > 0)])
    scPA = np.arctan2(negdRa,negdDec)*(u.rad).to("deg")
    scPA2= np.arctan2(posdRa,posdDec)*(u.rad).to("deg")

    if (scPA > 0) and (scPA < 180):
        myPA = scPA
    else:
        myPA = scPA2
        
#    import pdb;pdb.set_trace()

    #scanlen = np.max([nkotf.XSize,nkotf.YSize])/2.0
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
        print 'Span = ',span
        print nXpix,nYpix

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
                          avgRA,avgDEC,w,RRmap,nkotf)

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
        
#        if i % 100 == 0:
#            print np.sum(Coverage.time)/(len(Coverage.myFOVy)*i/fSamp)
#            pdb.set_trace()
            
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

def Single_Scan(skyobj, nkotf=None, date=None, doPlot=False, pwv=4,elMin=40,
                **kwargs):

    goodMin, mydate,goodStart = find_good_times(elMin=elMin,skyobj=skyobj,date=date)
    time_per_day    = (np.max(goodMin) - np.min(goodMin))*u.min
    if nkotf == None:
        nkotf = nkotf_parameters(**kwargs)
    TimeAr, ScanX, ScanY = nkotf_scan(XSize=nkotf.XSize,YSize=nkotf.YSize,PA=nkotf.PA,
                                      Tilt=nkotf.Tilt, Step=nkotf.Step, Speed=nkotf.Speed,
                                      fSamp=nkotf.fSamp)
    time_per_scan   = (np.max(TimeAr) - np.min(TimeAr))*u.s
    tpS_w_overhead  = (time_per_scan*1.1).to("min")

    ### If I wanted to try to simulate exactly the overhead (pointings, focus)
    ### I might want to take breaks every hour and every 2-3 hours.
    ### But I will hold off on this.
    ScansperHour    = int(((1.0*u.hour)/tpS_w_overhead).decompose().value)
    ScansperDay     = int(((time_per_day)/tpS_w_overhead).decompose().value)
    Coverage = None
    NEFD1mm, eta1mm, Tau1mm = values_by_band(band="1mm",pwv=pwv)
    NEFD2mm, eta2mm, Tau2mm = values_by_band(band="2mm",pwv=pwv)

    ScanNum = 0
    ### The starting Minute:
    MinStart = np.min(goodMin) + (ScanNum*tpS_w_overhead).to("min").value
    ### The starting time (date + hour,min,second):
    tStart = create_mytime(mydate,MinStart)
    tStart = goodStart     # This should be the better way to do it...
    if doPlot == True:
        plot_radecs(RaDecs)
    Coverage = coverage_map(tStart,nkotf,skyobj,TimeAr, ScanX, ScanY,
                            Coverage=Coverage,pwv=pwv)
    print "Scan "+str(ScanNum)+" of "+str(ScansperDay)

    Coverage.weight1mm *= eta1mm/(NEFD1mm**2)
    Coverage.weight2mm *= eta2mm/(NEFD2mm**2)
    nzwt1mm = (Coverage.weight1mm > 0); nzwt2mm = (Coverage.weight2mm > 0)
    Coverage.noise1mm[nzwt1mm] = Coverage.weight1mm[nzwt1mm]**(-0.5)
    Coverage.noise2mm[nzwt2mm] = Coverage.weight2mm[nzwt2mm]**(-0.5)

    return Coverage

def Multiple_Scans(tObs, skyobj, nkotf=None, date=None, doPlot=False, pwv=4,elMin=40,
                   **kwargs):

    goodMin, mydate,goodStart = find_good_times(elMin=elMin,skyobj=skyobj,date=date)
    time_per_day    = (np.max(goodMin) - np.min(goodMin))*u.min
    if nkotf == None:
        nkotf = nkotf_parameters(**kwargs)
    TimeAr, ScanX, ScanY = nkotf_scan(XSize=nkotf.XSize,YSize=nkotf.YSize,PA=nkotf.PA,
                                      Tilt=nkotf.Tilt, Step=nkotf.Step, Speed=nkotf.Speed,
                                      fSamp=nkotf.fSamp)
    time_per_scan   = (np.max(TimeAr) - np.min(TimeAr))*u.s
    tpS_w_overhead  = (time_per_scan*1.1).to("min")

    ### If I wanted to try to simulate exactly the overhead (pointings, focus)
    ### I might want to take breaks every hour and every 2-3 hours.
    ### But I will hold off on this.
    ScansperHour    = int(((1.0*u.hour)/tpS_w_overhead).decompose().value)
    ScansperDay     = int(((time_per_day)/tpS_w_overhead).decompose().value)
    Coverage = None
    NEFD1mm, eta1mm, Tau1mm = values_by_band(band="1mm",pwv=pwv)
    NEFD2mm, eta2mm, Tau2mm = values_by_band(band="2mm",pwv=pwv)

    for ScanNum in range(ScansperDay):
        ### The starting Minute:
        MinStart = np.min(goodMin) + (ScanNum*tpS_w_overhead).to("min").value
        ### The starting time (date + hour,min,second):
        tStart = create_mytime(mydate,MinStart)
        tStart = goodStart + (ScanNum*tpS_w_overhead) # date object
        if doPlot == True:
            plot_radecs(RaDecs)
        Coverage = coverage_map(tStart,nkotf,skyobj,TimeAr, ScanX, ScanY,
                                Coverage=Coverage,pwv=pwv)
        print "Scan "+str(ScanNum)+" of "+str(ScansperDay)

    Coverage.weight1mm *= eta1mm/(NEFD1mm**2)
    Coverage.weight2mm *= eta2mm/(NEFD2mm**2)
    nzwt1mm = (Coverage.weight1mm > 0); nzwt2mm = (Coverage.weight2mm > 0)
    Coverage.noise1mm[nzwt1mm] = Coverage.weight1mm[nzwt1mm]**(-0.5)
    Coverage.noise2mm[nzwt2mm] = Coverage.weight2mm[nzwt2mm]**(-0.5)

    return Coverage

def make_fits(Coverage,target="Object"):

    header = Coverage.w.to_header()

    hdu1 = fits.PrimaryHDU(Coverage.time,header=header)
    hdu1.header.append(("Title","Time Map"))
    hdu1.header.append(("Target",target))
    ### Weight Maps:
    hdu2 = fits.ImageHDU(Coverage.weight1mm)
    hdu2.header = header
    hdu2.name = 'Weight_Map_1mm'
    hdu2.header.append(("Title","Weight Map 1mm"))
    hdu2.header.append(("Target",target))
    hdu2.header.append(("XTENSION","Second"))
    hdu2.header.append(("SIMPLE","T")) 
    hdu2.verify('fix')
    hdu3 = fits.ImageHDU(Coverage.weight2mm)
    hdu3.header = header
    hdu3.name = 'Weight_Map_2mm'
    hdu3.header.append(("Title","Weight Map 1mm"))
    hdu3.header.append(("Target",target))
    hdu3.header.append(("XTENSION","Second"))
    hdu3.header.append(("SIMPLE","T")) 
    hdu3.verify('fix')
    ### Noise Maps:
    hdu4 = fits.ImageHDU(Coverage.noise1mm)
    hdu4.header = header
    hdu4.name = 'Noise_Map_1mm'
    hdu4.header.append(("Title","Noise Map 1mm"))
    hdu4.header.append(("Target",target))
    hdu4.header.append(("XTENSION","Second"))
    hdu4.header.append(("SIMPLE","T")) 
    hdu4.verify('fix')
    hdu5 = fits.ImageHDU(Coverage.noise2mm)
    hdu5.header = header
    hdu5.name = 'Noise_Map_2mm'
    hdu5.header.append(("Title","Noise Map 1mm"))
    hdu5.header.append(("Target",target))
    hdu5.header.append(("XTENSION","Second"))
    hdu5.header.append(("SIMPLE","T")) 
    hdu5.verify('fix')

    hdu1.header.add_history("Coverage maps made on "+Coverage.date_made+ ".")
    hdu1.header.add_history("Coverage maps are for observations on "+Coverage.date_obs+ ".")
    hdulist = fits.HDUList([hdu1,hdu2,hdu3,hdu4,hdu5])
    hdulist.info()
    filename="Coverage_Maps_"+target+".fits"
    fullpath = os.path.join(cwd,filename)
    hdulist.writeto(fullpath,overwrite=True,output_verify="exception")

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




#################################################################################




#################################################################################



#################################################################################
#################################################################################
##########################                             ##########################
#############                                                       #############
#############                   PLOTTING ROUTINES                   #############
#############                                                       #############
##########################                             ##########################
#################################################################################
#################################################################################



#################################################################################




#################################################################################




#################################################################################

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
#    AllAvgEC1 = 1.1
#    AllAvgEC2 = 1.1
    
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

#def add_target_text(Coverage):

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
#    plt.text(myxx,myyy,'hello')
#    print myxx,myyy
#    import pdb; pdb.set_trace()

def hist_pas(Coverage,dpi=200,filename="Position_Angle_Histogram_",
             addname="Object",myfontsize=15,format='png'):

    fig = plt.figure(dpi=200,figsize=(8,8))
    ax = fig.add_axes([0.10, 0.3, 0.85, 0.65])
    ax1 = fig.add_axes([0.10, 0.10, 0.85, 0.10])
    # N is the count in each bin, bins is the lower-limit of the bin
    N, bins, patches = ax.hist(Coverage.scanPA,bins='auto')
    binsize = np.median(bins - np.roll(bins,1))
    binEl  = np.array([]);  paStart= np.array([]); paStop = np.array([])
    
    for pa in bins:
        gi = (Coverage.scanPA >= pa) & (Coverage.scanPA < pa+binsize)
        myEls = Coverage.scanel[gi].value
        binEl  = np.append(binEl,np.mean(myEls))
        paStart = np.append(paStart,np.min(Coverage.scanstart[gi]).value)
        paStop  = np.append(paStop,np.max(Coverage.scanstop[gi]).value)
        #ncou   = len(myEls)
        
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
        mytime =  datetime.datetime.strptime(mystart, '%Y-%m-%d %H:%M:%S.%f')
        myutc  = mytime.strftime('%H:%M')+' UTC'
        ax.text(pa,yy/2.0+4,myutc,color='black',fontsize=myfontsize,rotation=-90)

    startdt = datetime.datetime.strptime(Coverage.scanstart[0].value, '%Y-%m-%d %H:%M:%S.%f')
    startday= startdt.strftime('%Y-%m-%d')
    myxloc = np.min(bins)*0.6 + np.max(bins)*0.4
    
    ax.text(myxloc,0.9*N, startday,color='black',fontsize=myfontsize)
    
    units="degrees"
    cb1 = cb.ColorbarBase(ax1, cmap=cmap,norm=norm,
                          orientation='horizontal')
    ax.set_title("Histogram of Position Angles on "+addname,fontsize=myfontsize)
    ax.set_xlabel("Position Angle (degrees)",fontsize=myfontsize)
    ax.set_ylabel("Number of Scans",fontsize=myfontsize)
    ax1.set_xlabel("Average Elevation",fontsize=myfontsize)
    fullbase = os.path.join(cwd,filename)
    fulleps = fullbase+addname+'.eps'; fullpng = fullbase+addname+'.png'
    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')
    plt.clf()

#    import pdb;pdb.set_trace()

#############################################################################
    
def ind_plots_cov(Coverage,map,filename="NIKA2_Coverage_map",target="Object",
                  myfontsize=15,mytitle="my map",addname="_quantity",
                  units='(units)',band="1mm",cblim=False,addtext=False,dpi=200,
                  secObj=None,thiObj=None,fouObj=None,format='png'):

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
    fullbase = os.path.join(cwd,filename)
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
                  secObj=None,thiObj=None,fouObj=None,format='png'):

    myfontsize=15
    ind_plots_cov(Coverage,Coverage.time,filename=filename,target=target,
                  mytitle="Time Map; ",addname="_time",units='seconds',
                  myfontsize=myfontsize,secObj=secObj,thiObj=thiObj,
                  fouObj=fouObj,format=format)
    
    ind_plots_cov(Coverage,Coverage.weight1mm,filename=filename,target=target,
                  mytitle="Weight Map, 1mm; ",addname="_weight1mm",units="mJy/beam $^{-2}$",
                  myfontsize=myfontsize,band="1mm",secObj=secObj,thiObj=thiObj,
                  fouObj=fouObj,format=format)

    ind_plots_cov(Coverage,Coverage.weight2mm,filename=filename,target=target,
                  mytitle="Weight Map, 2mm; ",addname="_weight2mm",units="mJy/beam $^{-2}$",
                  myfontsize=myfontsize,band="2mm",secObj=secObj,thiObj=thiObj,
                  fouObj=fouObj,format=format)

    ind_plots_cov(Coverage,Coverage.noise1mm,filename=filename,target=target,
                  mytitle="Noise Map, 1mm; ",addname="_noise1mm",units="mJy/beam",
                  myfontsize=myfontsize,band="1mm",cblim=True,addtext=True,secObj=secObj,thiObj=thiObj,
                  fouObj=fouObj,format=format)

    ind_plots_cov(Coverage,Coverage.noise2mm,filename=filename,target=target,
                  mytitle="Noise Map, 2mm; ",addname="_noise2mm",units="mJy/beam",
                  myfontsize=myfontsize,band="2mm",cblim=True,addtext=True,secObj=secObj,thiObj=thiObj,
                  fouObj=fouObj,format=format)

def plot_visibility(mydate,mysky,Coverage=None,mylabel="Target",dpi=200,elMin=40,
                    filename = "Visibility_Chart",utcoffset=2.0*u.hour,format='png'):

    bestdate = find_transit(mysky,mydate,utcoffset=utcoffset)

#    midnight = create_mytime(mydate,0,utcoffset=2.0*u.hour)
    mid_date = datetime.datetime.strptime(bestdate.value[0], "%Y-%m-%d %H:%M:%S.%f")
    npts = 1000
    date_arr = np.array([mid_date + datetime.timedelta(hours=i*24.0/npts - 12.0) for i in xrange(npts)])
    delta_bestdate = np.linspace(-12, 12, npts)*u.hour
#    delta_bestdate = np.linspace(-12, 12, npts)*u.hour
#    mytimes  = bestdate + delta_bestdate
    mytimes  = bestdate + delta_bestdate
    myframe  = apc.AltAz(obstime=mytimes, location=thirtym)
    sunaltazs = apc.get_sun(mytimes).transform_to(myframe)
    moonaltazs= apc.get_moon(mytimes).transform_to(myframe)
    objaltazs = mysky.transform_to(myframe)

    GoodEl = (objaltazs.alt.value > elMin)
    elStart= np.min(date_arr[GoodEl])
    elStop = np.max(date_arr[GoodEl])

    deltaT = 24.0*u.hour / npts
    bt30 = len(mytimes[(objaltazs.alt.value > 30.0)]) * deltaT.value
    bt40 = len(mytimes[(objaltazs.alt.value > 40.0)]) * deltaT.value
    bt50 = len(mytimes[(objaltazs.alt.value > 50.0)]) * deltaT.value
    bteM = len(mytimes[(objaltazs.alt.value > elMin)]) * deltaT.value

    plt.figure(1,dpi=dpi,figsize=(8,8));    plt.clf();    fig1,ax1 = plt.subplots()

#    mydates = matplotlib.dates.date2num(mytimes)

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
#    import pdb;pdb.set_trace()

    myylim = ax1.get_ylim(); ax1.set_ylim(myylim)
    ax1.plot(myxlim,[elMin,elMin],'--',color ='b',
             label="{:.1f}".format(bteM)+" hrs above "+str(int(elMin)))
    ax1.plot_date([elStart,elStart],[0,90],'--',color='b')
    ax1.plot_date([elStop,elStop],[0,90],'--',color='b')

    
    gi = np.array([])
    if Coverage != None:
        for start,stop in zip(Coverage.scanstart,Coverage.scanstop):
            mgi = np.where(((mytimes > start) & (mytimes < stop)) == True)
            gi = np.append(gi,mgi)

#        import pdb; pdb.set_trace()
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
    fullbase = os.path.join(cwd,filename)
    fulleps = fullbase+'.eps'; fullpng = fullbase+'.png'
    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')

def plot_radecs(RaDecs,span=600,format='png'):
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
    
    filename = "RaDec_map_v2";fullbase = os.path.join(cwd,filename)
    fulleps = fullbase+'.eps'; fullpng = fullbase+'.png'
    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')

def plot_altaz(ScanAz,ScanEl,TestRaDec,format='png'):

    plt.figure()
    plt.plot(ScanAz.value,ScanEl.value,'.')
    filename = "AltAz_map_v2";fullbase = os.path.join(cwd,filename)
    fulleps = fullbase+'.eps'; fullpng = fullbase+'.png'
    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')
        
    plt.figure()
    plt.plot(TestRaDec.ra.value,TestRaDec.dec.value,'.')
    filename = "TestRaDec_map_v2";fullbase = os.path.join(cwd,filename)
    fulleps = fullbase+'.eps'; fullpng = fullbase+'.png'
    if format == 'png':
        plt.savefig(fullpng,format='png')
    else:
        plt.savefig(fulleps,format='eps')
