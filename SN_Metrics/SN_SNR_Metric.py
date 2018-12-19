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
        
    def run(self, dataSlice, slicePoint=None):
        
        goodFilters = np.in1d(dataSlice['filter'],self.filterNames)
        dataSlice = dataSlice[goodFilters]
        if dataSlice.size == 0:
            return None
        dataSlice.sort(order=self.mjdCol)
        
        time = dataSlice[self.mjdCol]-dataSlice[self.mjdCol].min()
        r = []
        shift = 10 # SN DayMax: current date - shift days
        seasons = [self.season]
        self.max_rf_phase = -20
        self.min_rf_phase = 40
        
        self.margin = (1.+self.config['Observations']['z']) * (self.max_rf_phase-self.min_rf_phase) / 365.25
        
        if self.season == -1:
            seasons = np.unique(dataSlice[self.seasonCol])

        r = []
        """
        result_queue = multiprocessing.Queue()
        for season in seasons:
            p=multiprocessing.Process(name='Subprocess-'+str(season),target=self.Ana_Season,args=(dataSlice,season,shift,season,result_queue))
            p.start()
           

        resultdict = {}
        for season in seasons:
            resultdict.update(result_queue.get())
        
        for p in multiprocessing.active_children():
            p.join()

        for season in seasons:
             r.append(resultdict[season])  

        """
        for season in seasons:
            proc_season = self.Ana_Season(dataSlice,season,shift)
            if proc_season is not None:
                r.append(proc_season)
                
        return np.concatenate(r)

    def Ana_Season(self,dataSlice,season,shift,j=-1,output_q=None):
   
        idx = dataSlice[self.seasonCol] == season
        sel = dataSlice[idx]
        if len(sel) == 0:
            return None
        fieldRA = np.mean(sel[self.RaCol])
        fieldDec =  np.mean(sel[self.DecCol])
        Nvisits = np.median(sel[self.nexpCol]/2.) # one visit = 2 exposures
        m5 = np.mean(sel[self.m5Col])
        exptime = np.median(sel[self.exptimeCol])
        sel.sort(order = self.mjdCol)
        mjds = sel[self.mjdCol]
        mjd_min, mjd_max = np.min(mjds), np.max(mjds)
        season_length = mjd_max-mjd_min
        dates = np.arange(mjd_min, mjd_max+1., 1.)
        T0_lc = dates-shift
        cadence = np.mean(mjds[1:]-mjds[:-1])
        band = np.unique(sel[self.filterCol])[0]
        diff_time = dates[:,np.newaxis]-mjds
        time_for_lc = -T0_lc[:,None]+mjds
        flag = diff_time >=0
        
        print('T0',T0_lc)
        print('MJDS',mjds)
        print(time_for_lc)
        print('hello',sel[[self.m5Col,self.mjdCol]])
        m5_vals = np.tile(sel[self.m5Col],(len(time_for_lc),1))
        mjd_vals = np.tile(sel[self.mjdCol],(len(time_for_lc),1))
        
        res_new = None
        fluxes_tot= []
        for ib,interp in enumerate(self.names_ref):
            name = self.names_ref[ib]
            fluxes,snr = self.SNR_new(time_for_lc,m5_vals,flag,ib)
            fluxes_tot.append(fluxes)
            if res_new is None:
                res_new = np.array(snr, dtype=[('SNR_'+name,'f8')])
            else:
                res_new = rf.append_fields(res_new,'SNR_'+name,snr)

        res_new = rf.append_fields(res_new,'DayMax',T0_lc)
        res_new = rf.append_fields(res_new,'MJD',dates)
        res_new = rf.append_fields(res_new,'m5_eff',np.mean(np.ma.array(m5_vals,mask=~flag), axis=1))
        
        global_info = [(fieldRA,fieldDec,band,season_length,cadence,m5,mjd_min,season,Nvisits,exptime)]*len(res_new)
        names = ['fieldRA','fieldDec','band','season_length','Cadence','m5','MJD_min','season','Nvisits','ExposureTime']
        global_info = np.rec.fromrecords(global_info,names = names)
        res_new = rf.append_fields(res_new,names,[global_info[name] for name in names])
     
        self.Plot_new(fluxes_tot,mjd_vals,flag,res_new,self.config['names_ref'])
        """
        print(time_for_lc[:,-1])
        print(np.where(diff_time >=0.))
        """
        idx = diff_time >=0
      
        r=[]
        names = ['fieldRA','fieldDec','band','MJD_min','season_length','Nvisits','m5','ExposureTime','Cadence']
        ro = [fieldRA,fieldDec,band,mjd_min,mjd_max-mjd_min,Nvisits,m5,exptime,cadence]
        if self.display:
            import pylab as plt 
            plt.ion()
            fig, ax = plt.subplots(ncols=1, nrows=2)
            fig.canvas.draw()
            fig.suptitle('Ra = '+str(np.round(fieldRA,2))+' Dec = '+str(np.round(fieldDec,2))+' \n '+band+' band - season '+str(season))
            fieldid = str(np.round(fieldRA,2))+'_'+str(np.round(fieldDec,2))
        T0 = dates-shift
    
        for io,id in enumerate(idx):
            T0 = dates[io]-shift
            #if T0 > mjd_min+self.margin and len(sel[id])>=2:
            if T0 > -666.:
                sel_tot,fluxes_tot,res, names = self.Fill_Vals(sel[id],dates,io,T0,ro, names)
                if res is not None:
                    r.append(res)
        
                if self.display:
                    self.Plot(fig,ax,sel_tot,fluxes_tot,T0,shift,cadence,dates[io],dates[-1],
                              self.config['names_ref'],
                              np.rec.fromrecords(r,names = names),
                              fieldid,io)
            
        res = np.rec.fromrecords(r,names = names)
        if output_q is not None:
            output_q.put({j : res})
        else:
            print(res[['m5_eff','SNR_SNSim','SNR_SNCosmo']],res_new[['m5_eff','SNR_SNSim','SNR_SNCosmo']])
            print(len(res),len(res_new))
            print(res.dtype)
            print(res_new.dtype)
            return res

    def Fill_Vals(self, sel_ref,dates,io,T0,r_ref, names_ref):

        ro = r_ref.copy()
        names = names_ref.copy()
        dict_vals = dict(zip(names_ref,r_ref))
        #cadence_eff = np.mean(sel[self.mjdCol][id][1:]-sel[self.mjdCol][id][:-1])
        #m5_eff = np.mean(sel[self.m5Col][id])
        cadence_eff = np.mean(sel_ref[self.mjdCol][1:]-sel_ref[self.mjdCol][:-1])
        m5_eff = np.mean(sel_ref[self.m5Col])
        #this is to estimate the results with a "perfect" cadence

        #fakes = np.array(np.arange(np.min(sel_ref[self.mjdCol]),np.max(sel_ref[self.mjdCol]),cadence_eff),dtype=[(self.mjdCol,'f8')])
        fakes = np.array(np.arange(np.min(sel_ref[self.mjdCol]),np.max(sel_ref[self.mjdCol]),3.),dtype=[(self.mjdCol,'f8')])
        fakes = rf.append_fields(fakes,self.m5Col,[m5_eff]*len(fakes))
        
        ro+=[dates[io],T0,cadence_eff,m5_eff]
        names+=['MJD','DayMax','Cadence_eff','m5_eff']
                
        sel_tot = {}
        fluxes_tot = {}

        for season in np.unique(sel_ref['season']):
            idxa = sel_ref['season'] == season
            sel_season = sel_ref[idxa]
            names.append('season')
            ro.append(season)
            for ib,interp in enumerate(self.names_ref):
                    name = self.names_ref[ib]
                    sel,fluxes,snr = self.SNR(sel_season,T0,ib)
                    fluxes_tot[name] = fluxes
                    sel_tot[name] = sel
                    """
                    ro .append(band)
                    if 'band' not in names:
                    names.append('band')
                    """
                    ro.append(snr)
                    names.append('SNR_'+self.names_ref[ib])

            """
            sel,fluxes,snr = self.SNR(fakes,T0,ib)
            fluxes_tot[name+'_fake'] = fluxes
            sel_tot[name+'_fake'] = sel
            ro.append(snr)
            names.append('SNR_'+self.names_ref[ib]+'_fake')
            """
        return sel_tot,fluxes_tot,tuple(ro),names

    def SNR(self,sel,T0,ib):
        
        fluxes = self.lim_sn.fluxes[ib](np.copy(sel[self.mjdCol])-T0)
        flux_5sigma = self.lim_sn.mag_to_flux[ib](np.copy(sel[self.m5Col]))
        snr = np.sum(fluxes**2/flux_5sigma**2)
        snr = 5.*np.sqrt(snr)
                        
        return sel,fluxes,snr

    def SNR_new(self, time_lc, m5_vals, flag,ib):

        fluxes = self.lim_sn.fluxes[ib](time_lc)
        flux_5sigma = self.lim_sn.mag_to_flux[ib](m5_vals)
        snr = fluxes**2/flux_5sigma**2
        
        return fluxes,5.*np.sqrt(np.sum(snr*flag,axis =1))

    def Plot_new(self,fluxes, mjd,flag,snr, names):

        import pylab as plt
        plt.ion()
        fig, ax = plt.subplots(ncols=1, nrows=2)
        fig.canvas.draw()
      
        colors = ['b','r']
        myls = ['-','--']

        mjd_ma = np.ma.array(mjd,mask=~flag)
        #for ib in range(len(fluxes)):
        fluxes_ma = {}
        for ib, name in enumerate(names):
            fluxes_ma[name] = np.ma.array(fluxes[ib], mask=~flag)
        """
        ib=0
        flux = fluxes[ib]
        fluxes_ma = np.ma.array(flux, mask=~flag)
        """
        key = names[0]
        jmax = len(fluxes_ma[key])
        for j in range(len(fluxes_ma[key])):
            for ib, name in enumerate(names):
                ax[0].errorbar(mjd_ma[j],fluxes_ma[name][j],marker='s',color=colors[ib],ls = myls[ib])
                ax[1].plot(snr['MJD'][:j],snr['SNR_'+name][:j],color=colors[ib])
            fig.canvas.flush_events()
            if j != jmax-1:
                ax[0].clear()
        
                
    def Plot (self,fig,ax,sel_tot,fluxes_tot,T0,shift,cadence,current_mjd,mjd_max,name,snr,fieldid,inum):
        
        colors = ['b','r']
        myls = ['-','--']
        mfc = ['b','None']
        tot_label = []
        fontsize = 12
        fluxes = []
        for key, flux in fluxes_tot.items():
            fluxes.append(flux)
        min_fluxes = np.min(np.concatenate(fluxes))
        max_fluxes = np.max(np.concatenate(fluxes))
        
        tot_label.append(ax[0].errorbar([T0,T0],[min_fluxes,max_fluxes],color='k',label = 'DayMax'))
        tot_label.append(ax[0].errorbar([current_mjd,current_mjd],[min_fluxes,max_fluxes],color='k',ls='--', label ='Current MJD'))

        ib = -1
        for key, sel in sel_tot.items():
            if 'fake' not in key:
                ib+=1
                tot_label.append(ax[0].errorbar(sel[self.mjdCol],fluxes_tot[key],marker='s',color=colors[ib],ls = myls[ib],label = key))
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
                #ax[1].plot(snr['MJD'],snr['SNR_'+val+'_fake'],color=colors[i],ls = '--')
            else:
               ax[1].plot(snr['MJD'],snr['SNR_'+val],color=colors[i],label = val)
               #ax[1].plot(snr['MJD'],snr['SNR_'+val+'_fake'],color=colors[i],label = val+'_fake',ls='--')

        ax[1].set_xlabel('MJD',fontsize = fontsize)
        ax[1].set_ylabel('SNR',fontsize = fontsize)
        ax[1].legend()
        fig.canvas.flush_events()

        #fig.savefig('Test_figs/'+fieldid+'_'+str(inum)+'.png')

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
