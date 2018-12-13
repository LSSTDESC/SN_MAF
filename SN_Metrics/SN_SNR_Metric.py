import numpy as np
from lsst.sims.maf.metrics import BaseMetric
from Coadd_Stacker import CoaddStacker
import healpy as hp
import numpy.lib.recfunctions as rf

class SNMetric(BaseMetric):
    """
    Measure SN-SNR as a function of time.
    """
    def __init__(self, metricName='SNMetric',
                 mjdCol='observationStartMJD', RaCol='fieldRA',DecCol='fieldDec',
                 filterCol='filter', m5Col='fiveSigmaDepth',exptimeCol='visitExposureTime',
                 nightCol='night',obsidCol='observationId',nexpCol='numExposures',
                 vistimeCol='visitTime',coadd=True,lim_sn = None,names_ref = None,
                 uniqueBlocks=False, config=None,**kwargs):
    
        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.RaCol= RaCol
        self.DecCol = DecCol
        self.exptimeCol = exptimeCol
        self.seasonCol = 'season'
        self.nightCol = nightCol
        self.obsidCol = obsidCol
        self.nexpCol = nexpCol
        self.vistimeCol = vistimeCol
        
        cols = [self.nightCol,self.m5Col,self.filterCol,self.mjdCol,self.obsidCol,self.nexpCol,self.vistimeCol, self.exptimeCol,self.seasonCol]
        if coadd:
            cols+=['coadd']
        super(SNMetric, self).__init__(col=cols,metricDtype = 'object',metricName=metricName, **kwargs)
        
        self.filterNames = np.array(['u','g','r','i','z','y'])
        self.config = config

        # this is for output
        """
        save_status = config['Output']['save']
        outdir = config['Output']['directory']
        prodid = config['ProductionID']
        """
        # sn parameters
        sn_parameters = config['SN parameters']
        
        self.field_type = config['Observations']['fieldtype']
        self.season = config['Observations']['season']
        #self.season = season
        area = 9.6 #survey_area in sqdeg - 9.6 by default for DD
        if  self.field_type == 'WFD':
            # in that case the survey area is the healpix area
            area = hp.nside2pixarea(config['Pixelisation']['nside'],degrees=True)
        
        # Load the reference Li file

        #self.Li = np.load(config['Reference File'])
        self.lim_sn = lim_sn
        self.names_ref = names_ref

    def run(self, dataSlice, slicePoint=None):
        # Cut down to only include filters in correct wave range.
        
        #print("attention",type(dataSlice),dataSlice.dtype)
        goodFilters = np.in1d(dataSlice['filter'],self.filterNames)
        dataSlice = dataSlice[goodFilters]
        if dataSlice.size == 0:
            return None
        dataSlice.sort(order=self.mjdCol)
        #print('dataslice',np.unique(dataSlice[['fieldRA','fieldDec','season','filter']]),dataSlice.dtype)
        time = dataSlice[self.mjdCol]-dataSlice[self.mjdCol].min()
        r = []
        shift = 10 # SN DayMax: current date - shift days
        seasons = [self.season]
       
        if self.season == -1:
            seasons = np.unique(dataSlice[self.seasonCol])
        for season in seasons:
            r.append(self.Ana_Season(dataSlice,season,shift))

        """
        import pylab as plt
        res = np.concatenate(r)
        for season in np.unique(res['season']):
            ii = res['season']==season
            sel = res[ii]
            plt.plot(sel['MJD'],sel['SNR_SNSim'],color = 'b')
            plt.plot(sel['MJD'],sel['SNR_SNCosmo'],color = 'r')
        plt.show()
        """                         
        return np.concatenate(r)

    def Ana_Season(self,dataSlice,season,shift):
        names = ['fieldRA','fieldDec','season','band']
        idx = dataSlice[self.seasonCol] == season
        sel = dataSlice[idx]
        fieldRA = np.mean(sel[self.RaCol])
        fieldDec =  np.mean(sel[self.DecCol])
        band  = np.unique(sel[self.filterCol])[0]
        
        sel.sort(order = self.mjdCol)
        mjd_min = np.min(sel[self.mjdCol])
        mjd_max = np.max(sel[self.mjdCol])
        dates = np.arange(mjd_min, mjd_max+1., 1.)
        diff_time = dates[:,np.newaxis]-sel[self.mjdCol]
        cadence = np.mean(sel[self.mjdCol][1:]-sel[self.mjdCol][:-1])
        idx = diff_time >=0
        r=[]
        for io,id in enumerate(idx):
            T0 = dates[io]-shift
            ro = []
            if T0 > mjd_min:
                ro += [fieldRA,fieldDec,season,band]
                ro+=[dates[io],T0,cadence]
                if 'MJD' not in names:
                    names+=['MJD','DayMax','Cadence']
                for ib,interp in enumerate(self.names_ref):
                    fluxes = self.lim_sn.fluxes[ib](T0-sel[self.mjdCol][id])
                    flux_5sigma = self.lim_sn.mag_to_flux[ib](sel[self.m5Col][id])
                    #print(fluxes,flux_5sigma)
                    snr = np.sum(fluxes**2/flux_5sigma**2)
                    snr = 5.*np.sqrt(snr)
                        
                    ro.append(snr)
                    if 'SNR_'+self.names_ref[ib] not in names:
                        names.append('SNR_'+self.names_ref[ib])
            if ro:
                r.append(tuple(ro))
       
        res = np.rec.fromrecords(r,names = names)
        return res
