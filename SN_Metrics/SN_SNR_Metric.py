import numpy as np
from lsst.sims.maf.metrics import BaseMetric
from Coadd_Stacker import CoaddStacker
import healpy as hp
import numpy.lib.recfunctions as rf
import multiprocessing

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
        
    def run(self, dataSlice, seasons_=None,slicePoint=None):
        
        goodFilters = np.in1d(dataSlice['filter'],self.filterNames)
        dataSlice = dataSlice[goodFilters]
        if dataSlice.size == 0:
            return None
        dataSlice.sort(order=self.mjdCol)
        
        time = dataSlice[self.mjdCol]-dataSlice[self.mjdCol].min()
        #r = []
        shift = 10 # SN DayMax: current date - shift days
        
        seasons = self.season
        if seasons_ is not None:
            seasons=seasons_
            
        if self.season == -1:
            seasons = np.unique(dataSlice[self.seasonCol])

        proc_season = self.Ana_Season(dataSlice,seasons,shift)
        """
        if proc_season is not None:
            r.append(proc_season)
        """
        return proc_season
    
    def Ana_Season(self,dataSlice,seasons,shift,j=-1,output_q=None):

        sel = None
        """
        periods = {}
        cadence = {}
        season_length = {}
        mjd_min = {}
        """
        rv = []
        for season in seasons:
            idx =(dataSlice[self.seasonCol] == season)
            slice_sel = dataSlice[idx]
            slice_sel.sort(order=self.mjdCol)
            mjds_season = slice_sel[self.mjdCol]
            #periods[season]= (np.min(mjds_season),np.max(mjds_season))
            cadence  = np.mean(mjds_season[1:]-mjds_season[:-1])
            mjd_min = np.min(mjds_season)
            mjd_max = np.max(mjds_season)
            season_length = mjd_max-mjd_min
            rv.append((float(season),cadence,season_length,mjd_min,mjd_max))
            if sel is None:
                sel = dataSlice[idx]
            else:
                sel = np.concatenate((sel, dataSlice[idx]))

        self.info_season = np.rec.fromrecords(rv, names=['season','cadence','season_length','MJD_min','MJD_max'])
        
        if len(sel) == 0:
            return None
        fieldRA = np.mean(sel[self.RaCol])
        fieldDec =  np.mean(sel[self.DecCol])
        Nvisits = np.median(sel[self.nexpCol]/2.) # one visit = 2 exposures
        m5 = np.mean(sel[self.m5Col])
        exptime = np.median(sel[self.exptimeCol])
        sel.sort(order = self.mjdCol)
        mjds = sel[self.mjdCol]
        dates = None
        #for key, val in periods.items():
        for val in self.info_season:
            if dates is None:
                dates = np.arange(val['MJD_min'], val['MJD_max']+1., 1.)
            else:
                dates = np.concatenate((dates, np.arange(val['MJD_min'], val['MJD_max']+1., 1.)))
        
        T0_lc = dates-shift
        
        band = np.unique(sel[self.filterCol])[0]
        diff_time = dates[:,np.newaxis]-mjds
        time_for_lc = -T0_lc[:,None]+mjds
        flag = (diff_time >=0)&(diff_time<=200.)
        
        m5_vals = np.tile(sel[self.m5Col],(len(time_for_lc),1))
        mjd_vals = np.tile(sel[self.mjdCol],(len(time_for_lc),1))
        season_vals = np.tile(sel[self.seasonCol],(len(time_for_lc),1))
    
        fluxes_tot, snr= self.SNR(time_for_lc,m5_vals,flag,season_vals)
        
        _,idx = np.unique(snr['season'],return_inverse=True)
        infos = self.info_season[idx]
        vars_info = ['cadence','season_length','MJD_min']
        snr = rf.append_fields(snr,vars_info,[infos[name] for name in vars_info])
        snr = rf.append_fields(snr,'DayMax',T0_lc)
        snr = rf.append_fields(snr,'MJD',dates)
        snr = rf.append_fields(snr,'m5_eff',np.mean(np.ma.array(m5_vals,mask=~flag), axis=1))
        global_info = [(fieldRA,fieldDec,band,m5,Nvisits,exptime)]*len(snr)
        names = ['fieldRA','fieldDec','band','m5','Nvisits','ExposureTime']
        global_info = np.rec.fromrecords(global_info,names = names)
        snr = rf.append_fields(snr,names,[global_info[name] for name in names])
        if self.display:
            self.Plot(fluxes_tot,mjd_vals,flag,snr,T0_lc,dates)
            
        if output_q is not None:
            output_q.put({j : snr})
        else:
            return snr

    def SNR(self, time_lc, m5_vals, flag,season_vals):

        seasons = np.ma.array(season_vals,mask=~flag)
        fluxes_tot = {}
        snr_tab = None
        
        for ib,name in enumerate(self.names_ref):
            fluxes = self.lim_sn.fluxes[ib](time_lc)
            if name not in fluxes_tot.keys():
                fluxes_tot[name] = fluxes
            else:
                fluxes_tot[name] = np.concatenate((fluxes_tot[name],fluxes))
                
            flux_5sigma = self.lim_sn.mag_to_flux[ib](m5_vals)
            snr = fluxes**2/flux_5sigma**2
            snr_season = 5.*np.sqrt(np.sum(snr*flag,axis =1))
            if snr_tab is None:
                snr_tab = np.asarray(np.copy(snr_season),dtype=[('SNR_'+name,'f8')])
            else:
                snr_tab= rf.append_fields(snr_tab,'SNR_'+name,np.copy(snr_season))
        snr_tab = rf.append_fields(snr_tab,'season',np.mean(seasons, axis = 1))
        
        return fluxes_tot,snr_tab
    
    def Plot(self,fluxes, mjd,flag,snr,T0_lc, dates):

        import pylab as plt
        plt.ion()
        fig, ax = plt.subplots(ncols=1, nrows=2)
        fig.canvas.draw()
      
        colors = ['b','r']
        myls = ['-','--']
        mfc = ['b','None']
        tot_label = []
        fontsize = 12
        mjd_ma = np.ma.array(mjd,mask=~flag)
        fluxes_ma = {}
        for key, val in fluxes.items():
            fluxes_ma[key] = np.ma.array(val, mask=~flag)

        key = list(fluxes.keys())[0]
        jmax = len(fluxes_ma[key])
        tot_label=[]
        tot_label_snr = []
        min_flux = []
        max_flux =  []
        for j in range(jmax):
             
            for ib,name in enumerate(fluxes_ma.keys()):
                tot_label.append(ax[0].errorbar(mjd_ma[j],fluxes_ma[name][j],marker='s',color=colors[ib],ls = myls[ib],label=name))
                #tot_label_snr.append(ax[1].errorbar(snr['MJD'][:j],snr['SNR_'+name][:j],color=colors[ib],label=name,ls='None',marker='.',ms=3.))
                tot_label_snr.append(ax[1].errorbar(snr['MJD'][:j],snr['SNR_'+name][:j],color=colors[ib],label=name))
                min_flux.append(np.min(fluxes_ma[name][j]))
                max_flux.append(np.max(fluxes_ma[name][j]))
            min_fluxes = np.min(min_flux)
            max_fluxes = np.max(max_flux)
            tot_label.append(ax[0].errorbar([T0_lc[j],T0_lc[j]],[min_fluxes,max_fluxes],color='k',label = 'DayMax'))
            tot_label.append(ax[0].errorbar([dates[j],dates[j]],[min_fluxes,max_fluxes],color='k',ls='--', label ='Current MJD'))
            fig.canvas.flush_events()
            if j != jmax-1:
                ax[0].clear()
                tot_label=[]
                tot_label_snr = []

        labs = [l.get_label() for l in tot_label]
        ax[0].legend(tot_label, labs, ncol=1,loc='best',prop={'size':fontsize},frameon=False)
        ax[0].set_ylabel('Flux [e.sec$^{-1}$]',fontsize = fontsize)

        ax[1].set_xlabel('MJD',fontsize = fontsize)
        ax[1].set_ylabel('SNR',fontsize = fontsize)
        ax[1].legend()
        labs = [l.get_label() for l in tot_label_snr]
        ax[1].legend(tot_label_snr, labs, ncol=1,loc='best',prop={'size':fontsize},frameon=False)
        for i in range(2):
            ax[i].tick_params(axis='x', labelsize=fontsize)
            ax[i].tick_params(axis='y', labelsize=fontsize)
