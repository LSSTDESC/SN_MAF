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
from sn_cadence_tools import Reference_Data, Generate_Fake_Observations
from scipy import interpolate

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

    for band in bands:
        sql_i = sqlconstraint+' AND '
        sql_i += 'filter = "%s"' % (band)

        lim_sn[band] = Reference_Data(
            config['Li file'], config['Mag_to_flux file'], band, z)

        metric[band] = module.SNMetric(config=config, coadd=config['Observations']
                                       ['coadd'], lim_sn=lim_sn[band], names_ref=config['names_ref'])
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
    res = None
    for band, val in bdict.items():
        print('ici band', band)
        res_b = np.concatenate(val.metricValues[~val.metricValues.mask])
        res_b = np.unique(res_b)
        if res is None:
            res = res_b
        else:
            res = np.concatenate((res, res_b))

    print('bands', np.unique(res['band']))
    """
    for (Ra, Dec, season) in np.unique(res[['fieldRA', 'fieldDec', 'season']]):
        idx = (res['fieldRA'] == Ra) & (
            res['fieldDec'] == Dec) & (res['season'] == season)
        sel = res[idx]
        Plot_SNR(Ra, Dec, season, sel, config, metric, z)
    """

    frac_obs = Fraction_Observation(res, config, metric)
    print(frac_obs)
    # mbg.writeAll()
    # mbg.plotAll(closefigs=False)
    # mbg.plot()
    plt.show()


def Plot_SNR(Ra, Dec, season, data, config, metric, z, draw_fakes=True):
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
        ifig = int(ib/2)
        jfig = int(ib % 2)

        ax = axa[ifig][jfig]
        if draw_fakes:  # superimpose fake obs with a "perfect" cadence
            fake_obs = Get_Fake_Obs(sel, config, band)  # simulate fakes
            # run the same metric as simulated observations
            resb = metric[band].run(
                fake_obs[fake_obs['filter'] == band], seasons_=[1])
        sel.sort(order='MJD')

        # Draw results
        for io, sim in enumerate(config['names_ref']):
            tot_label.append(ax.errorbar(
                sel['MJD'], sel['SNR_'+sim], ls='-', color=colors[io], label=sim))
            if draw_fakes:
                tot_label.append(ax.errorbar(
                    resb['MJD'], resb['SNR_'+sim], ls='--', color=colors[io], label=sim+'_fake'))

        if ifig == ncols-1:
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


def Fraction_Observation(data, config, metric):
    """
    Estimate the fraction of time (during a season)
    a supernovae can be observed ie with
    a minimal SNR (for each band)
    """
    r = []
    for (ra, dec, season) in np.unique(data[['fieldRA', 'fieldDec', 'season']]):
        idx = (data['fieldRA'] == ra) & (
            data['fieldDec'] == dec) & (data['season'] == season)
        sela = data[idx]
        # produce fake observations
        fake_obs = Get_Fake_Obs(sel, config, band)  # simulate fakes
        # run the same metric as simulated observations
        resb = metric[band].run(
            fake_obs[fake_obs['filter'] == band], seasons_=[1])
        resb.sort(order='MJD')
        sel.sort(order='MJD')
        r = [ra, dec, season, band]
        names = ['Ra', 'Dec', 'season', 'band']
        for sim in config['names_ref']:
            fakes = interpolate.interp1d(resb['MJD'], resb['SNR_'+sim])
            obs = interpolate.interp1d(sel['MJD'], sel['SNR_'+sim])
            mjd_min = np.max([np.min(sel['MJD']), np.min(resb['MJD'])])
            mjd_max = np.min([np.max(sel['MJD']), np.max(resb['MJD'])])
            mjd = np.arange(mjd_min, mjd_max, 1.)
            print(fakes(mjd))
            print(obs(mjd))
            diff_res = obs(mjd)-fakes(mjd)
            idx = diff_res >= 0
            print(len(diff_res[idx]), len(diff_res[idx])/len(diff_res))
            r += [len(diff_res[idx])/len(diff_res)]
            names += ['frac_obs_'+sim]

        return np.rec.fromrecords(r, names=names)


def Get_Fake_Obs(sel, config, band):
    """
    Simulation of fake observations
    Needed information (cadence, m5, mjd, season length, Nvisits ...) are
    extracted from simulations
    """

    config_fake = yaml.load(open(config['Fake_file']))
    # m5_ref = dict(zip('grizy', [23.27, 24.58, 24.22, 23.65, 22.78, 22.0]))
    cadence = np.mean(sel['cadence'])
    mjd_min = np.mean(sel['MJD_min'])
    season_length = np.mean(sel['season_length'])
    Nvisits = np.median(sel['Nvisits'])
    m5 = np.median(sel['m5'])
    Tvisit = 30.
    # cadence = 3.

    config_fake['bands'] = [band]
    config_fake['Cadence'] = [cadence]
    config_fake['MJD_min'] = mjd_min
    config_fake['season_length'] = season_length
    config_fake['Nvisits'] = [Nvisits]
    m5_nocoadd = m5-1.25*np.log10(float(Nvisits)*Tvisit/30.)
    config_fake['m5'] = [m5_nocoadd]
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
