import numpy as np
from lsst.sims.maf.metrics import BaseMetric
from SN_Simulation import SN_Simulation
from Coadd_Stacker import CoaddStacker

class SNMetric(BaseMetric):
    """
    Measure how many time series meet a given time and filter distribution requirement.
    """
    def __init__(self, metricName='SNMetric',
                 mjdCol='observationStartMJD', filterCol='filter', m5Col='fiveSigmaDepth',
                 units='', redshift=0.,
                 Tmin = -20., Tmax = 60., Nbetween=7, Nfilt={'u':0,'g':10,'r':10,'i':10,'z':5,'y':5}, Nfilt_obs=5,Tless = -5., Nless=1,
                 Tmore = 30., Nmore=1, peakGap=2., snrCut=10., singleDepthLimit=23.,
                 resolution=5., badval=-666,
                 uniqueBlocks=False, config=None,**kwargs):
        """
        redshift = redshift of the SN.  Used to scale observing dates to SN restframe.
        Tmin = the minimum day to consider the SN.
        Tmax = the maximum to consider.
        Nbetween = the number of observations to demand between Tmin and Tmax
        Nfilt = number of unique filters that must observe the SN above the snrCut
        Tless = minimum time to consider 'near peak'
        Tmore = max time to consider 'near peak'
        Nless = number of observations to demand before Tless
        Nmore = number of observations to demand after Tmore
        peakGap = maximum gap alowed between observations in the 'near peak' time
        snrCut = require snr above this limit when counting Nfilt XXX-not yet implemented
        singleDepthLimit = require observations in Nfilt different filters to be this
        deep near the peak.  This is a rough approximation for the Science Book
        requirements for a SNR cut.  Ideally, one would import a time-variable SN SED,
        redshift it, and make filter-keyed dictionary of interpolation objects so the
        magnitude of the SN could be calculated at each observation and then use the m5col
        to compute a SNR.
        resolution = time step (days) to consider when calculating observing windows
        uniqueBlocks = should the code count the number of unique sequences that meet
        the requirements (True), or should all sequences that meet the conditions
        be counted (False).

        The filter centers are shifted to the SN restframe and only observations
        with filters between 300 < lam_rest < 900 nm are included

        In the science book, the metric demands Nfilt observations above a SNR cut.
        Here, we demand Nfilt observations near the peak with a given singleDepthLimt.
        """
        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.dateCol = 'expDate'
        self.fieldRA='fieldRA'
        self.fieldDec='fieldDec'
        self.ditheredRA='ditheredRA'
        self.ditheredDec='ditheredDec'
        self.visitTime='visitExpTime'
        self.finSeeing='finSeeing'
        self.rawSeeing='rawSeeing'
        self.moonPhase='moonPhase'
        self.airmass='airmass'
        self.filtSkyBrightness='filtSkyBrightness'

        """
        super(SNMetric, self).__init__(col=[self.mjdCol, self.m5Col, self.filterCol, self.dateCol,self.fieldRA,self.fieldDec, self.ditheredRA,self.ditheredDec,self.visitTime,self.finSeeing,self.rawSeeing,self.moonPhase,self.airmass,self.filtSkyBrightness],
                                              metricName=metricName, units=units, badval=badval,
                                              **kwargs)
        """
        super(SNMetric, self).__init__(col=['night','fiveSigmaDepth','filter',self.mjdCol,'observationId','airmass','numExposures', 'visitTime', 'visitExposureTime','season','coadd'],metricName=metricName, units=units, badval=badval,**kwargs)
        self.redshift = redshift
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Nbetween = Nbetween
        self.Nfilt = Nfilt
        self.Nfilt_obs = Nfilt_obs
        self.Tless = Tless
        self.Nless = Nless
        self.Tmore = Tmore
        self.Nmore = Nmore
        self.peakGap = peakGap
        self.snrCut = snrCut
        self.resolution = resolution
        self.uniqueBlocks = uniqueBlocks
        self.filterNames = np.array(['u','g','r','i','z','y'])
        self.config = config
        print('hello', self.config)

        # load cosmology
        cosmo_par = config['Cosmology']
        # load telescope
        tel_par = config['Instrument']

        # this is for output

        save_status = config['Output']['save']
        outdir = config['Output']['directory']
        prodid = config['ProductionID']
        # sn parameters
        sn_parameters = config['SN parameters']

        simu_config = config['Simulator']
        display_lc = config['Display']

        names = dict(zip(['band', 'mjd', 'rawSeeing', 'sky', 'exptime',
                      'moonPhase', 'Ra', 'Dec', 'Nexp', 'fiveSigmaDepth',
                      'seeing', 'airmass', 'night', 'season', 'pixarea',
                      'pixRa', 'pixDec'],
                     ['band', 'mjd', 'seeingFwhm500', 'sky', 'exptime',
                      'moonPhase', 'Ra', 'Dec', 'numExposures',
                      'fiveSigmaDepth', 'seeingFwhmEff', 'airmass',
                      'night', 'season', 'pixarea',
                     'pixRa', 'pixDec']))

        self.simu = SN_Simulation(cosmo_par, tel_par, sn_parameters,
                             save_status, outdir, prodid,
                             simu_config, display_lc, names=names,nproc=config['Multiprocessing']['nproc'])

        

    def run(self, dataSlice, slicePoint=None):
        # Cut down to only include filters in correct wave range.
        
        #print("attention",type(dataSlice),dataSlice.dtype)
        goodFilters = np.in1d(dataSlice['filter'],self.filterNames)
        dataSlice = dataSlice[goodFilters]
        if dataSlice.size == 0:
            return (self.badval, self.badval,self.badval)
        dataSlice.sort(order=self.mjdCol)
        print('dataslice',np.unique(dataSlice[['fieldRA','fieldDec','season']]),dataSlice.dtype)
        time = dataSlice[self.mjdCol]-dataSlice[self.mjdCol].min()

        #print 'time in Rest frame',len(time),time[191],time[192],time[1380],time[1381],time[1382]
        # Creat time steps to evaluate at
        finetime = np.arange(0.,np.ceil(np.max(time)),self.resolution)
        #print 'finetime',finetime
        #index for each time point
        ind = np.arange(finetime.size)
        #index for each time point + Tmax - Tmin
        right = np.searchsorted( time, finetime+self.Tmax-self.Tmin, side='right')
        left = np.searchsorted(time,finetime, side='left')
        #print 'left and right',left,right
        # Demand enough visits in window  
        good = np.where( (right - left) > self.Nbetween)[0]
        #print 'right minus left',right,left,right - left,'good',good
        ind = ind[good]
        right = right[good]
        left = left[good]
        #print 'index',ind
        result = 0
        # Record the maximum gap near the peak (in rest-frame days)
        maxGap = []
        # Record the total number of observations in a sequence.
        Nobs = []
        # Record the number of observations betwenn Tless and Tmore
        Nobs_peak = []
        outputToASCII = True
        
 
        maxGap = np.array(maxGap)
        Nobs=np.array(Nobs)
        Nobs_peak=np.array(Nobs_peak)
        #print 'hello pal',Nobs
        return {'result':result, 'maxGap':maxGap, 'Nobs':Nobs, 'Nobs_peak':Nobs_peak}

    def reduceMedianMaxGap(self,data):
        """The median maximum gap near the peak of the light curve """
        result = np.median(data['maxGap'])
        if np.isnan(result):
            result = self.badval
        return result
    def reduceNsequences(self,data):
        """The number of sequences that met the requirements """
        return data['result']
    def reduceMedianNobs(self,data):
        """Median number of observations covering the entire light curve """
        result = np.median(data['Nobs'])
        if np.isnan(result):
            result = self.badval
        return result

    def reduceMedianNobs_peak(self,data):
        """Median number of observations covering the entire light curve """
        result = np.median(data['Nobs_peak'])
        if np.isnan(result):
            result = self.badval
        return result
