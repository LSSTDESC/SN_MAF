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

        self.display = config['Display_Processing']
        
    def run(self, dataSlice, slicePoint=None):
        
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
        self.max_rf_phase = -20
        self.min_rf_phase = 40
        #print(self.config)
        self.margin = (1.+self.config['Observations']['z']) * (self.max_rf_phase-self.min_rf_phase) / 365.25
        if self.season == -1:
            seasons = np.unique(dataSlice[self.seasonCol])
        for season in seasons:
            proc_season = self.Ana_Season(dataSlice,season,shift)
            if proc_season is not None:
                r.append(proc_season)
                
        return np.concatenate(r)

    def Ana_Season(self,dataSlice,season,shift):
   
        idx = dataSlice[self.seasonCol] == season
        sel = dataSlice[idx]
        if len(sel) == 0:
            return None
        fieldRA = np.mean(sel[self.RaCol])
        fieldDec =  np.mean(sel[self.DecCol])
        band  = np.unique(sel[self.filterCol])[0]
        Nvisits = np.median(sel[self.nexpCol]/2.) # one visit = 2 exposures
        m5 = np.mean(sel[self.m5Col])
        exptime = np.median(sel[self.exptimeCol])
        sel.sort(order = self.mjdCol)
        mjd_min = np.min(sel[self.mjdCol])
        mjd_max = np.max(sel[self.mjdCol])
        dates = np.arange(mjd_min, mjd_max+1., 1.)
        diff_time = dates[:,np.newaxis]-sel[self.mjdCol]
        cadence = np.mean(sel[self.mjdCol][1:]-sel[self.mjdCol][:-1])
        names = ['fieldRA','fieldDec','season','band','MJD_min','season_length','Nvisits','m5','ExposureTime','Cadence']
        idx = diff_time >=0
        r=[]
        if self.display:
            import pylab as plt 
            plt.ion()
            fig, ax = plt.subplots(ncols=1, nrows=2)
            fig.canvas.draw()
            fig.suptitle('Ra = '+str(fieldRA)+' Dec = '+str(fieldDec)+' season = '+str(season))
            fieldid = str(np.round(fieldRA,2))+'_'+str(np.round(fieldDec,2))
            
        for io,id in enumerate(idx):
            T0 = dates[io]-shift
            #print('go',T0,dates[io],dates[-1])
            ro = []
            if T0 > mjd_min+self.margin and len(sel[id])>=2:
                cadence_eff = np.mean(sel[self.mjdCol][id][1:]-sel[self.mjdCol][id][:-1])
                m5_eff = np.mean(sel[self.m5Col][id])
                ro += [fieldRA,fieldDec,season,band,mjd_min,mjd_max-mjd_min,Nvisits,m5,exptime,cadence]
                ro+=[dates[io],T0,cadence_eff,m5_eff]
                if 'MJD' not in names:
                    names+=['MJD','DayMax','Cadence_eff','m5_eff']
                
                sel_tot = []
                fluxes_tot = []
                for ib,interp in enumerate(self.names_ref):
                    fluxes = self.lim_sn.fluxes[ib](T0-sel[self.mjdCol][id])
                    fluxes_tot.append(fluxes)
                    sel_tot.append(sel[id])
                    
                    flux_5sigma = self.lim_sn.mag_to_flux[ib](sel[self.m5Col][id])
                    snr = np.sum(fluxes**2/flux_5sigma**2)
                    snr = 5.*np.sqrt(snr)
                        
                    ro.append(snr)
                    if 'SNR_'+self.names_ref[ib] not in names:
                        names.append('SNR_'+self.names_ref[ib])
                if ro:
                    r.append(ro)
                if self.display:
                    self.Plot(fig,ax,sel_tot,fluxes_tot,T0,shift,cadence,dates[io],dates[-1],
                              self.config['names_ref'],
                              np.rec.fromrecords(r,names = names),
                              fieldid,io)
            
        #print(r)
        res = np.rec.fromrecords(r,names = names)
        return res

    def Plot (self,fig,ax,sel_tot,fluxes_tot,T0,shift,cadence,current_mjd,mjd_max,name,snr,fieldid,inum):
        #import pylab as plt
        colors = ['b','r']
        myls = ['-','--']
        mfc = ['b','None']
        tot_label = []
        fontsize = 12
        min_fluxes = np.min(np.concatenate(fluxes_tot))
        max_fluxes = np.max(np.concatenate(fluxes_tot))
        
        tot_label.append(ax[0].errorbar([T0,T0],[min_fluxes,max_fluxes],color='k',label = 'DayMax'))
        tot_label.append(ax[0].errorbar([current_mjd,current_mjd],[min_fluxes,max_fluxes],color='k',ls='--', label ='Current MJD'))

        for ib,sel in enumerate(sel_tot):
            fluxes = fluxes_tot[ib]
            tot_label.append(ax[0].errorbar(sel[self.mjdCol],fluxes,marker='s',color=colors[ib],ls = myls[ib],label = name[ib]))
            #test = np.arange(np.min(sel[self.mjdCol]),current_mjd,cadence)
            #tot_label.append(ax[0].errorbar(test, self.lim_sn.fluxes[ib](T0-test),color=colors[ib],marker='o',label = name[ib]+' regular cadence',ls='None',markerfacecolor=mfc[ib]))
            
        fig.canvas.flush_events()
        labs = [l.get_label() for l in tot_label]
        ax[0].legend(tot_label, labs, ncol=1,loc='best',prop={'size':fontsize},frameon=False)
        #ax[0].set_xlabel('MJD',fontsize = fontsize)
        ax[0].set_ylabel('Flux [e.sec$^{-1}$]',fontsize = fontsize)
       

        for i, val in enumerate(name):
            if mjd_max - current_mjd != 0.:
                ax[1].plot(snr['MJD'],snr['SNR_'+val],color=colors[i])
            else:
               ax[1].plot(snr['MJD'],snr['SNR_'+val],color=colors[i],label = val) 

        ax[1].set_xlabel('MJD',fontsize = fontsize)
        ax[1].set_ylabel('SNR',fontsize = fontsize)
        ax[1].legend()
        fig.canvas.flush_events()

        fig.savefig('Test_figs/'+fieldid+'_'+str(inum)+'.png')

        if mjd_max - current_mjd != 0.:
            ax[0].clear()
            tot_label= []
        
        """
        plt.pause(1)
        plt.close()
        
        plt.draw()
        plt.pause(1)
        plt.close()
        """
