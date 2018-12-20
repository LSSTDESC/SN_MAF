import matplotlib.pyplot as plt
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.db as db
import lsst.sims.maf.utils as utils
import argparse
import time
import yaml
from importlib import import_module
import sqlite3
import numpy as np
from SN_Cadence_Tools import Reference_Data,Generate_Fake_Observations

parser = argparse.ArgumentParser(
    description='Run a SN metric from a configuration file')
parser.add_argument('config_filename',
                    help='Configuration file in YAML format.')

def run(config_filename):
    # YAML input file.
    config = yaml.load(open(config_filename))
    #print(config)
    outDir ='Test' # this is for MAF

    # grab the db filename from yaml input file 
    dbFile=config['Observations']['filename']

    """
    conn = sqlite3.connect(dbFile)
    cur = conn.cursor()
    table_name='Proposal'
    result = cur.execute("PRAGMA table_info('%s')" % table_name).fetchall()
    print('Results',result)

    cur.execute("SELECT * FROM Proposal")
    rows = cur.fetchall()
    for row in rows:
        print(row)
    print('end')
    cur.execute('PRAGMA TABLE_INFO({})'.format('ObsHistory'))
    
    names = [tup[1] for tup in cur.fetchall()]
    print(names)
    """
    opsimdb = db.OpsimDatabase(dbFile)
    version = opsimdb.opsimVersion
    propinfo, proptags = opsimdb.fetchPropInfo()
    print('proptags and propinfo',proptags,propinfo)

    # grab the fieldtype (DD or WFD) from yaml input file
    fieldtype = config['Observations']['fieldtype']
    
    module = import_module(config['Metric'])

    slicer = slicers.HealpixSlicer(nside=config['Pixelisation']['nside'])
    
    sqlconstraint = opsimdb.createSQLWhere(fieldtype, proptags)

    bundles = []
    names = []
    lim_sn = {}
    bands  = config['Observations']['bands']
    z = config['Observations']['z']
    metric ={}
    for band in bands:
        sql_i = sqlconstraint+' AND '
        sql_i += 'filter = "%s"' % (band)
        
        lim_sn[band]= Reference_Data(config['Li file'],config['Mag_to_flux file'],band,z)
        
        metric[band]=module.SNMetric(config=config,coadd=config['Observations']['coadd'], lim_sn = lim_sn[band],names_ref=config['names_ref'])
        bundles.append(metricBundles.MetricBundle(metric[band], slicer, sql_i))
        names.append(band)

    bdict = dict(zip(names,bundles))
   
    resultsDb = db.ResultsDb(outDir='None')
    mbg =  metricBundles.MetricBundleGroup(bdict, opsimdb,
                                      outDir=outDir, resultsDb=resultsDb)
    
    result = mbg.runAll()  
    
    for band, val in bdict.items():
        res = np.concatenate(val.metricValues[~val.metricValues.mask])
        res = np.unique(res)
        for (Ra,Dec,season) in np.unique(res[['fieldRA','fieldDec','season']]):
            idx = (res['fieldRA'] == Ra)&(res['fieldDec'] == Dec)&(res['season'] == season)
            sel = res[idx]
            Plot_SNR(Ra,Dec,season,band,sel,config,metric)
    
        
    # Let us display the results
    """
    for band, val in bdict.items():
        lim_sn[band].Plot_Cadence_Metric(val.metricValues[~val.metricValues.mask])
        lim_sn[band].Plot_Hist_zlim(config['names_ref'],val.metricValues[~val.metricValues.mask])
    """
    #mbg.writeAll()
    #mbg.plotAll(closefigs=False)
    #mbg.plot()
    plt.show()

def Plot_SNR(Ra,Dec,season,band,sel,config,metric,draw_fakes = True):
      colors = ['b','r']
      fontsize = 15

      tot_label=[]
      fig, ax = plt.subplots(ncols=1, nrows=1)
      fig.suptitle('Ra = '+str(np.round(Ra,2))+' Dec = '+str(np.round(Dec,2))+' \n '+band+' band - season '+str(int(season)),fontsize = fontsize)

      if draw_fakes:
          fake_obs = Get_Fake_Obs(sel,config,band)
          resb = metric[band].run(fake_obs[fake_obs['filter']==band],seasons_=[1])
      sel.sort(order='MJD')
      
      for io, sim in enumerate(config['names_ref']):
          tot_label.append(ax.errorbar(sel['MJD'],sel['SNR_'+sim],ls='-',color=colors[io], label = sim))
          if draw_fakes:
              tot_label.append(ax.errorbar(resb['MJD'],resb['SNR_'+sim],ls='--',color = colors[io],label = sim+'_fake'))
      ax.set_xlabel('MJD [day]',fontsize = fontsize)
      ax.set_ylabel('Signal-To-Noise ratio',fontsize = fontsize)
      ax.tick_params(axis='x', labelsize=fontsize)
      ax.tick_params(axis='y', labelsize=fontsize)
      labs = [l.get_label() for l in tot_label]
      ax.legend(tot_label, labs, ncol=1,loc='best',prop={'size':fontsize},frameon=False)

                           
def Get_Fake_Obs(sel,config,band):

    config_fake = yaml.load(open(config['Fake_file']))
    
    m5_ref = dict(zip('grizy',[23.27, 24.58, 24.22, 23.65, 22.78, 22.0]))

    cadence = np.mean(sel['cadence'])
    mjd_min = np.mean(sel['MJD_min'])
    season_length = np.mean(sel['season_length'])
    Nvisits = np.median(sel['Nvisits'])
    m5 = np.median(sel['m5'])
    Tvisit = 30.
    #cadence = 3.
    
    config_fake['bands'] = [band]
    config_fake['Cadence']= [cadence]
    config_fake['MJD_min'] = mjd_min
    config_fake['season_length'] = season_length
    config_fake['Nvisits'] = [Nvisits]
    m5_nocoadd = m5-1.25*np.log10(float(Nvisits)*Tvisit/30.)
    config_fake['m5']=[m5_nocoadd]
    fake_obs = Generate_Fake_Observations(config_fake).Observations
    return fake_obs

def main(args):
    print('running')
    time_ref = time.time()
    run(args.config_filename)
    print('Time', time.time()-time_ref)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
