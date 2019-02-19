import matplotlib.pyplot as plt
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.db as db
import argparse
import time
import yaml
from importlib import import_module
# import sqlite3
import numpy as np
from sn_maf.sn_tools.sn_cadence_tools import Reference_Data
import healpy as hp
import numpy.lib.recfunctions as rf

parser = argparse.ArgumentParser(
    description='Run a SN metric from a configuration file')
parser.add_argument('config_filename',
                    help='Configuration file in YAML format.')


def run(config_filename):
    # YAML input file.
    config = yaml.load(open(config_filename))
    # print(config)
    outDir = 'Test'  # this is for MAF

    # grab the db filename from yaml input file
    dbFile = config['Observations']['filename']

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
    # version = opsimdb.opsimVersion
    propinfo, proptags = opsimdb.fetchPropInfo()
    print('proptags and propinfo', proptags, propinfo)

    # grab the fieldtype (DD or WFD) from yaml input file
    fieldtype = config['Observations']['fieldtype']

    module = import_module(config['Metric'])

    slicer = slicers.HealpixSlicer(nside=config['Pixelisation']['nside'])

    sqlconstraint = opsimdb.createSQLWhere(fieldtype, proptags)

    bundles = []
    names = []
    lim_sn = {}
    bands = config['Observations']['bands']
    z = config['Observations']['z']
    metric = {}
    # processing. Band after band

    Ra_ref = 0.000
    Dec_ref = -2.308039

    for band in bands:
        sql_i = sqlconstraint+' AND '
        sql_i += 'filter = "%s"' % (band)
        #sql_i += ' AND abs(fieldRA-(%f))< %f' % (Ra_ref, 1.e-2)+' AND '
        #sql_i += 'abs(fieldDec-(%f))< %f' % (Dec_ref, 1.e-2)

        lim_sn[band] = Reference_Data(
            config['Li file'], config['Mag_to_flux file'], band, z)

        metric[band] = module.SNMetric(config=config, coadd=config['Observations']
                                       ['coadd'], lim_sn=lim_sn[band], names_ref=config['names_ref'], z=z)
        bundles.append(metricBundles.MetricBundle(metric[band], slicer, sql_i))
        names.append(band)

    bdict = dict(zip(names, bundles))

    resultsDb = db.ResultsDb(outDir='None')
    mbg = metricBundles.MetricBundleGroup(bdict, opsimdb,
                                          outDir=outDir, resultsDb=resultsDb)

    # result = mbg.runAll()
    mbg.runAll()
    # Plot the results:
    # SNR vs MJD for a SN with T0=MJD-10
    # Fake observations corresponding to a "perfect" cadence
    # can be superimposed

    # concatenate the results estimated per band
    metricValues = {}
    data_str = ['snr_obs', 'snr_fakes', 'detec_frac']

    for dstr in data_str:
        metricValues[dstr] = None

    for band, val in bdict.items():
        data = val.metricValues[~val.metricValues.mask]
        res = {}
        for dstr in data_str:
            res[dstr] = None
        for val in data:
            for dstr in data_str:
                if res[dstr] is None:
                    res[dstr] = val[dstr]
                else:
                    res[dstr] = np.concatenate((res[dstr], val[dstr]))

        for dstr in data_str:
            res[dstr] = np.unique(res[dstr])
            if metricValues[dstr] is None:
                metricValues[dstr] = res[dstr]
            else:
                metricValues[dstr] = np.concatenate(
                    (metricValues[dstr], res[dstr]))

    snr_obs = metricValues['snr_obs']
    snr_fakes = metricValues['snr_fakes']
    detec_frac = metricValues['detec_frac']

    for inum, (Ra, Dec, season) in enumerate(np.unique(snr_obs[['fieldRA', 'fieldDec', 'season']])):
        idx = (snr_obs['fieldRA'] == Ra) & (
            snr_obs['fieldDec'] == Dec) & (snr_obs['season'] == season)
        sel_obs = snr_obs[idx]
        idxb = (np.abs(snr_fakes['fieldRA'] - Ra) < 1.e-5) & (np.abs(
            snr_fakes['fieldDec'] - Dec) < 1.e-5) & (snr_fakes['season'] == season)
        sel_fakes = snr_fakes[idxb]
        SNRPlot(Ra, Dec, season, sel_obs, sel_fakes, config, metric, z)
        if inum >= 10:
            break

    print(detec_frac.dtype)

    DetecFracPlot(detec_frac, config['Pixelisation']
                  ['nside'], config['names_ref'])

    # frac_obs = Fraction_Observation(res, config, metric)
    # print(frac_obs)
    # mbg.writeAll()
    # mbg.plotAll(closefigs=False)
    # mbg.plot()
    plt.show()


def SNRPlot(Ra, Dec, season, data, data_fakes, config, metric, z, draw_fakes=True):
    """
    Signal-to-Ratio vs MJD plot
    SNR of  a SN with T0=MJD-10 days
    (x1,color) chosen in the input yaml file
    Fake observations can be superimposed
    One plot per field, per season.
    """

    colors = ['b', 'r']
    fontsize = 15
    bands_ref = 'ugrizy'
    id_band = [0, 1, 2, 3, 4, 5]
    bands_id = dict(zip(bands_ref, id_band))
    id_bands = dict(zip(id_band, bands_ref))
    bands = np.unique(data['band'])
    lista = sorted([bands_id[b] for b in bands])
    bands = [id_bands[jo] for jo in lista]
    n_bands = len(bands)
    # estimate the number of rows and columns depending on the number of bands
    ncols = 1
    nrows = 1
    if n_bands >= 2:
        ncols = 2
        nrows = int(n_bands/2+(n_bands % 2))

    figa, axa = plt.subplots(ncols=ncols, nrows=nrows, figsize=(15, 10))

    figa.suptitle('Ra = '+str(np.round(Ra, 2))+' Dec = '+str(np.round(Dec, 2)) +
                  ' \n '+' Season '+str(int(season))+' - z = '+str(z), fontsize=fontsize)
    for ib, band in enumerate(bands):
        tot_label = []
        idb = data['band'] == band
        sel = data[idb]
        idb = data_fakes['band'] == band
        sel_fakes = data_fakes[idb]
        sel.sort(order='MJD')
        sel_fakes.sort(order='MJD')
        ifig = int(ib/2)
        jfig = int(ib % 2)

        if nrows > 1:
            ax = axa[ifig][jfig]
        else:
            if ncols > 1:
                ax = axa[jfig]
            else:
                ax = axa

        # Draw results
        for io, sim in enumerate(config['names_ref']):
            tot_label.append(ax.errorbar(
                sel['MJD'], sel['SNR_'+sim], ls='-', color=colors[io], label=sim))
            if draw_fakes:
                tot_label.append(ax.errorbar(
                    sel_fakes['MJD'], sel_fakes['SNR_'+sim], ls='--', color=colors[io], label=sim+'_fake'))

        if ifig == nrows-1:
            ax.set_xlabel('MJD [day]', fontsize=fontsize)
        if jfig == 0:
            ax.set_ylabel('Signal-To-Noise ratio', fontsize=fontsize)

        ax.tick_params(axis='x', labelsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize)

        if ifig == 0 and jfig == 0:
            labs = [l.get_label() for l in tot_label]
            ax.legend(tot_label, labs, ncol=1, loc='best',
                      prop={'size': fontsize}, frameon=False)

        ax.text(0.9, 0.9, band, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes,
                fontsize=fontsize)


def DetecFracPlot(data, nside, names_ref):

    data_heal = GetHealpix(data, nside)
    npix = hp.nside2npix(nside)
    # print(data_heal)
    for band, season in np.unique(data_heal[['band', 'season']]):
        idx = (data_heal['band'] == band) & (data_heal['season'] == season)
        sel = data_heal[idx]
        for sim in names_ref:
            fig, ax = plt.subplots()
            hpxmap = np.zeros(npix, dtype=np.float)
            hpxmap[sel['healpixID']] += sel['frac_obs_'+sim]
            cmap = plt.cm.jet
            cmap.set_under('w')
            # remove max=200 and norm='hist' to get the DDFs
            median_value = np.median(sel['frac_obs_'+sim])
            plt.axes(ax)
            hp.mollview(hpxmap, min=0, max=1., cmap=cmap,
                        title='{} - season {} \n median: {}'.format(band, int(season), np.round(median_value, 2)), hold=True)

    plt.show()


def GetHealpix(data, nside):

    res = data.copy()
    npix = hp.nside2npix(nside)
    table = hp.ang2vec(data['fieldRA'], data['fieldDec'], lonlat=True)
    healpixs = hp.vec2pix(nside, table[:, 0], table[:, 1], table[:, 2])
    res = rf.append_fields(res, 'healpixID', healpixs)
    return res


def main(args):
    print('running')
    time_ref = time.time()
    run(args.config_filename)
    print('Time', time.time()-time_ref)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
