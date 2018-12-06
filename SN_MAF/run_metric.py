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
from scipy import interpolate

parser = argparse.ArgumentParser(
    description='Run a SN metric from a configuration file')
parser.add_argument('config_filename',
                    help='Configuration file in YAML format.')

def Lims(band, tab,SNR):

    lims = {}
    print(tab.dtype)
    for z in np.unique(tab['z']):
        lims[z] = {}
        idx = (tab['z'] == z)&(tab['band']=='LSST::'+band)
        sel = tab[idx]
        if len(sel) > 0:
            Li2 = np.sqrt(np.sum(sel['flux_e']**2))                                                               
            lim = 5. * Li2 / SNR
            lims[z][band]  = lim

    print('lims',lims)
    return lims

def Mesh(band, mag_to_flux,dt_range, mag_range,):
    dt = np.linspace(dt_range[0], dt_range[1], 50)
    m5 = np.linspace(mag_range[0], mag_range[1], 50)
    ida = mag_to_flux['band'] == band
    fa = interpolate.interp1d(mag_to_flux[ida]['m5'],mag_to_flux[ida]['flux'])
    f5 = fa(m5)
    F5,DT = np.meshgrid(f5, dt)
    M5,DT = np.meshgrid(m5, dt)
    metric = np.sqrt(DT) * F5

    return M5, DT, metric

def Plot(band,metricValues, Li, mag_to_flux,
         mag_range=(23., 27.5),dt_range=(0.5, 15.),
         target={# 'g': (26.91, 3.), # was 25.37                                            
                           'r': (26.5, 3.), # was 26.43                                               
                           'i': (26.16, 3.), # was 25.37      # could be 25.3 (400-s)                 
                           'z': (25.56, 3.), # was 24.68      # could be 25.1 (1000-s)                
                           'y': (24.68, 3.) },# was 24.72
         SNR=dict(zip([b for b in 'grizy'],                                             
                                [25., 25., 60., 35., 20.]))) :

    lims = []
    mags_to_fluxes = []
    M5_all = []
    DT_all = []
    metric_all = []
    
    for val in Li:
        lims.append(Lims(band, val, SNR[band]))
        
    for val in mag_to_flux:
        M5, DT, metric = Mesh(band, val,dt_range,mag_range)
        M5_all.append(M5)
        DT_all.append(DT)
        metric_all.append(metric)

    sorted_keys= []
    for i in range(len(lims)):
        sorted_keys.append(np.sort([k  for  k in  lims[i].keys()])[::-1])
        
    plt.figure()                                                                                               
    plt.imshow(metric, extent=(mag_range[0], mag_range[1],dt_range[0], dt_range[1]), aspect='auto', alpha=0.25)
    print(type(metricValues), len(np.copy(metricValues)),len(metricValues[~metricValues.mask]))
    for vval in metricValues:
            plt.plot(vval['m5_mean'],vval['cadence_mean'],'r+',alpha = 0.1)

    color = ['k','b']
    for kk,lim in enumerate(lims):
        fmt = {}                                                                                              
        ll = [lim[zz][band] for zz in sorted_keys[kk]]
        cs = plt.contour(M5_all[kk], DT_all[kk], metric_all[kk], ll, colors=color[kk])                                                       
        #dict_target_snsim=Get_target(cs,sorted_keys,cadence_ref,m5_exp)
        print('cs',cs.collections)
        for io,col in enumerate(cs.collections):
            print('col',sorted_keys[kk][io],col.get_segments())
        strs = ['$z=%3.1f$' % zz for zz in sorted_keys[kk]]
        for l,s in zip(cs.levels, strs):
            fmt[l] = s
        plt.clabel(cs, inline=True, fmt=fmt, fontsize=16, use_clabeltext=True)

        
    t = target.get(band, None)
    print(target, t)
    if t is not None:
        plt.plot(t[0], t[1],
                color='r', marker='*',
                markersize=15)
    plt.xlabel('$m_{5\sigma}$', fontsize=18)
    plt.ylabel(r'Observer frame cadence $^{-1}$ [days]', fontsize=18)
    plt.title('$%s$' % band.split(':')[-1], fontsize=18)
    plt.xlim(mag_range)
    plt.ylim(dt_range)
    plt.grid(1)     


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

    # print('slicer',slicer.pixArea,slicer.slicePoints['ra'])
    
    sqlconstraint = opsimdb.createSQLWhere(fieldtype, proptags)

    bundles = []
    #for band in 'grizy':
    names = []
    bands = 'r'
    for band in bands:
        sql_i = sqlconstraint+' AND '
        sql_i += 'filter = "%s"' % (band)
        #sql_i += ' AND '
        #sql_i +=  'season= "%s"' % (season)
        metric=module.SNMetric(config=config,coadd=config['Observations']['coadd'])
        bundles.append(metricBundles.MetricBundle(metric, slicer, sql_i))
        names.append(band)
        print('sql',sql_i)

    print('hello',len(bundles))
    #bdict = metricBundles.makeBundlesDictFromList(bundles)
    bdict = dict(zip(names,bundles))
    print('ahah',bdict.keys())
    """
    mb = metricBundles.MetricBundle(metric, slicer, sqlconstraint)
    
    mbD = {0:mb}

    resultsDb = db.ResultsDb(outDir='None')
    
    mbg =  metricBundles.MetricBundleGroup(mbD, opsimdb,
                                      outDir=outDir, resultsDb=resultsDb)

    """
    resultsDb = db.ResultsDb(outDir='None')
    mbg =  metricBundles.MetricBundleGroup(bdict, opsimdb,
                                      outDir=outDir, resultsDb=resultsDb)
    
    result = mbg.runAll()

    #for the plot, two files to be loaded
    Li = []
    mag_to_flux = []
    for val in config['Li file']:
        Li.append(np.load(val))
    for val in config['Mag_to_flux file']:
        mag_to_flux.append(np.load(val))
    SNR = dict(zip(bands,[25]))
    for band, val in bdict.items():
        Plot(band,val.metricValues[~val.metricValues.mask],Li,mag_to_flux)
                
        """
        for vval in val.metricValues:
            print(key, vval)
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
