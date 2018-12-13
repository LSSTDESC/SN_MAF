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
#from scipy import interpolate
#import numpy.lib.recfunctions as rf
#from scipy import interpolate
from SN_Cadence_Tools import Reference_Data
from SN_Utils import Generate_Fake_Observations

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
        #sql_i += ' AND '
        #sql_i +=  'season= "%s"' % (season)
        lim_sn[band]= Reference_Data(config['Li file'],config['Mag_to_flux file'],band,z)
        
        metric[band]=module.SNMetric(config=config,coadd=config['Observations']['coadd'], lim_sn = lim_sn[band],names_ref=config['names_ref'])
        bundles.append(metricBundles.MetricBundle(metric[band], slicer, sql_i))
        names.append(band)
       
        print('sql',sql_i)
    
    print('hello',len(bundles))
    bdict = dict(zip(names,bundles))
   
    resultsDb = db.ResultsDb(outDir='None')
    mbg =  metricBundles.MetricBundleGroup(bdict, opsimdb,
                                      outDir=outDir, resultsDb=resultsDb)
    
    result = mbg.runAll()

    config_fake = yaml.load(open(config['Fake_file']))
    fake_obs = Generate_Fake_Observations(config_fake).Observations
    print(config_fake)
    print(fake_obs.dtype)
    for band, val in bdict.items():
        res = np.concatenate(val.metricValues[~val.metricValues.mask])
        for (Ra,Dec,season) in np.unique(res[['fieldRA','fieldDec','season']]):
            idx = (res['fieldRA'] == Ra)&(res['fieldDec'] == Dec)&(res['season'] == season)
            sel = res[idx]
            cadence = np.mean(sel['Cadence'])
            config_fake['Cadence'][config_fake['bands'].index(band)] = cadence
            fake_obs = Generate_Fake_Observations(config_fake).Observations
            resb = metric[band].run(fake_obs[fake_obs['filter']==band])
            print(sel)
            print(resb)
    
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
    
def main(args):
    print('running')
    time_ref = time.time()
    run(args.config_filename)
    print('Time', time.time()-time_ref)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
