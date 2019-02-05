import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy.lib.recfunctions as rf
import matplotlib._cntr as cntr
import h5py
from astropy.table import Table, Column, vstack
from scipy.interpolate import griddata, interp2d, CloughTocher2DInterpolator
from sn_utils.utils.sn_telescope import Telescope


class Lims:
    def __init__(self, Li_files, mag_to_flux_files, band, SNR,
                 mag_range=(23., 27.5), dt_range=(0.5, 25.)):

        self.band = band
        self.SNR = SNR
        self.lims = []
        self.mag_to_flux = []
        self.mag_range = mag_range
        self.dt_range = dt_range

        for val in Li_files:
            self.lims.append(self.Get_Lims(self.band, np.load(val), SNR))
        for val in mag_to_flux_files:
            self.mag_to_flux.append(np.load(val))
        self.Interp()

    def Get_Lims(self, band, tab, SNR):

        lims = {}
        # print(tab.dtype)
        for z in np.unique(tab['z']):
            # lims[z] = {}
            idx = (tab['z'] == z) & (tab['band'] == 'LSST::'+band)
            idx &= (tab['flux_e'] > 0.)
            sel = tab[idx]

            if len(sel) > 0:
                Li2 = np.sqrt(np.sum(sel['flux_e']**2))
                lim = 5. * Li2 / SNR
                if z not in lims.keys():
                    lims[z] = {}
                lims[z][band] = lim

        # print('lims',lims)
        return lims

    def Mesh(self, mag_to_flux):
        dt = np.linspace(self.dt_range[0], self.dt_range[1], 100)
        m5 = np.linspace(self.mag_range[0], self.mag_range[1], 50)
        ida = mag_to_flux['band'] == self.band
        fa = interpolate.interp1d(
            mag_to_flux[ida]['m5'], mag_to_flux[ida]['flux_e'])
        f5 = fa(m5)
        F5, DT = np.meshgrid(f5, dt)
        M5, DT = np.meshgrid(m5, dt)
        metric = np.sqrt(DT) * F5

        return M5, DT, metric

    def Interp(self):
        M5_all = []
        DT_all = []
        metric_all = []

        for val in self.mag_to_flux:
            M5, DT, metric = self.Mesh(val)
            M5_all.append(M5)
            DT_all.append(DT)
            metric_all.append(metric)

        sorted_keys = []
        for i in range(len(self.lims)):
            sorted_keys.append(np.sort([k for k in self.lims[i].keys()])[::-1])
        figa, axa = plt.subplots()
        self.Points_Ref = []
        for kk, lim in enumerate(self.lims):
            fmt = {}
            ll = [lim[zz][self.band] for zz in sorted_keys[kk]]
            cs = axa.contour(M5_all[kk], DT_all[kk], metric_all[kk], ll)
            # print('cs',cs.collections)
            points_values = None
            for io, col in enumerate(cs.collections):
                if col.get_segments():
                    # print(sorted_keys[kk][io],'segment',col.get_segments())
                    myarray = col.get_segments()[0]
                    res = np.array(myarray[:, 0], dtype=[('m5', 'f8')])
                    res = rf.append_fields(res, 'cadence', myarray[:, 1])
                    res = rf.append_fields(
                        res, 'z', [sorted_keys[kk][io]]*len(res))
                    if points_values is None:
                        points_values = res
                    else:
                        points_values = np.concatenate((points_values, res))
            self.Points_Ref.append(points_values)
            """
            print('finally',restot)
            f = interpolate.interp2d(
                restot['m5'], restot['cadence'], restot['z'], kind='linear')
            self.Interpolate.append(f)
            """
        plt.close(figa)  # do not display

    def Interp_griddata(self, index, data):

        ref_points = self.Points_Ref[index]
        res = interpolate.griddata((ref_points['m5'], ref_points['cadence']), ref_points['z'], (
            data['m5_mean'], data['cadence_mean']), method='cubic')
        return res

    def Plot_Cadence_Metric(self, metricValues,
                            target={  # 'g': (26.91, 3.), # was 25.37
                                'r': (26.5, 3.),  # was 26.43
                                # was 25.37      # could be 25.3 (400-s)
                                'i': (26.16, 3.),
                                # was 24.68      # could be 25.1 (1000-s)
                                'z': (25.56, 3.),
                                'y': (24.68, 3.)}):  # was 24.72

        M5_all = []
        DT_all = []
        metric_all = []

        for val in self.mag_to_flux:
            M5, DT, metric = self.Mesh(val)
            M5_all.append(M5)
            DT_all.append(DT)
            metric_all.append(metric)

        sorted_keys = []
        for i in range(len(self.lims)):
            sorted_keys.append(np.sort([k for k in self.lims[i].keys()])[::-1])

        plt.figure()
        plt.imshow(metric, extent=(
            self.mag_range[0], self.mag_range[1], self.dt_range[0], self.dt_range[1]), aspect='auto', alpha=0.25)
        print(type(metricValues), len(np.copy(metricValues)),
              len(metricValues[~metricValues.mask]))
        restot = np.concatenate(metricValues)
        restot = np.unique(restot)
        # for vval in metricValues:
        plt.plot(restot['m5_mean'], restot['cadence_mean'], 'r+', alpha=0.9)

        color = ['k', 'b']
        for kk, lim in enumerate(self.lims):
            fmt = {}
            ll = [lim[zz][self.band] for zz in sorted_keys[kk]]
            cs = plt.contour(M5_all[kk], DT_all[kk],
                             metric_all[kk], ll, colors=color[kk])
            strs = ['$z=%3.1f$' % zz for zz in sorted_keys[kk]]
            for l, s in zip(cs.levels, strs):
                fmt[l] = s
            plt.clabel(cs, inline=True, fmt=fmt,
                       fontsize=16, use_clabeltext=True)

        t = target.get(self.band, None)
        if t is not None:
            plt.plot(t[0], t[1],
                     color='r', marker='*',
                     markersize=15)
        plt.xlabel('$m_{5\sigma}$', fontsize=18)
        plt.ylabel(r'Observer frame cadence $^{-1}$ [days]', fontsize=18)
        plt.title('$%s$' % self.band.split(':')[-1], fontsize=18)
        plt.xlim(self.mag_range)
        plt.ylim(self.dt_range)
        plt.grid(1)

    def Plot_Hist_zlim(self, names_ref, metricValues):
        r = []
        fontsize = 15
        restot = None
        colors = dict(zip(range(0, 4), ['r', 'k', 'b', 'g']))
        """
        for vval in metricValues:
            if restot is None:
                restot = np.copy(vval)
            else:
                    restot = np.concatenate((restot,np.copy(vval)))
         """
        restot = np.concatenate(metricValues)
        restot = np.unique(restot)
        fig, ax = plt.subplots(figsize=(8, 6))
        fig.suptitle(self.band + ' band', fontsize=fontsize)
        label = []
        xminv = []
        xmaxv = []
        for j, name in enumerate(names_ref):
            xminv.append(np.min(restot['zlim_'+name]))
            xmaxv.append(np.max(restot['zlim_'+name]))

        xmin = np.min(xminv)
        xmax = np.max(xmaxv)
        xstep = 0.025
        bins = np.arange(xmin, xmax+xstep, xstep)

        # print(restot[['band','m5_mean','cadence_mean','zlim_SNSim']])
        for j, name in enumerate(names_ref):
            label.append(
                name + '  $z_{med}$ = ' + str(np.median(np.round(restot['zlim_'+name], 2))))
            ax.hist(restot['zlim_'+name], range=[xmin, xmax],
                    bins=bins, histtype='step', color=colors[j], linewidth=2)

        ax.set_xlabel('$z_{lim}$', fontsize=fontsize)
        ax.set_ylabel(r'Number of Entries', fontsize=fontsize)
        # ax.set_xticks(np.arange(0.5,0.75,0.05))
        ax.tick_params(labelsize=fontsize)
        ax.grid()
        plt.legend(label, fontsize=fontsize-2., loc='upper left')
        # plt.grid(1)


class Reference_Data:
    def __init__(self, Li_files, mag_to_flux_files, band, z):

        self.band = band
        self.z = z
        self.fluxes = []
        self.mag_to_flux = []

        for val in Li_files:
            self.fluxes.append(self.Get_Interp_Fluxes(
                self.band, np.load(val), self.z))
        for val in mag_to_flux_files:
            self.mag_to_flux.append(
                self.Get_Interp_mag(self.band, np.load(val)))

    def Get_Interp_Fluxes(self, band, tab, z):

        lims = {}
        idx = (np.abs(tab['z'] - z) < 1.e-5) & (tab['band'] == 'LSST::'+band)
        sel = tab[idx]
        selc = np.copy(sel)
        difftime = (sel['time']-sel['DayMax'])
        selc = rf.append_fields(selc, 'deltaT', difftime)
        return interpolate.interp1d(selc['deltaT'], selc['flux_e'], bounds_error=False, fill_value=0.)

    def Get_Interp_mag(self, band, tab):
        print(tab.dtype)
        idx = tab['band'] == band
        sel = tab[idx]
        return interpolate.interp1d(sel['m5'], sel['flux_e'], bounds_error=False, fill_value=0.)


class Generate_Fake_Observations:
    """ Class to generate Fake observations
    Input
    ---------
    parameter file (filter, cadence, m5,Nseasons, ...)
    Production: Observations
    ---------
    recordarray of observations:
    MJD, Ra, Dec, band,m5,Nexp, ExpTime, Season
    """

    def __init__(self, config,
                 mjdCol='observationStartMJD', RaCol='fieldRA',
                 DecCol='fieldDec', filterCol='filter', m5Col='fiveSigmaDepth',
                 exptimeCol='visitExposureTime', nexpCol='numExposures', seasonCol='season'):

        # config = yaml.load(open(config_filename))
        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.RaCol = RaCol
        self.DecCol = DecCol
        self.exptimeCol = exptimeCol
        self.seasonCol = seasonCol
        self.nexpCol = nexpCol

        # now make fake obs
        self.make_fake(config)

    def make_fake(self, config):

        bands = config['bands']
        cadence = dict(zip(bands, config['Cadence']))
        shift_days = dict(
            zip(bands, [config['shift_days']*io for io in range(len(bands))]))
        m5 = dict(zip(bands, config['m5']))
        Nvisits = dict(zip(bands, config['Nvisits']))
        Exposure_Time = dict(zip(bands, config['Exposure_Time']))
        inter_season_gap = 300.

        Ra = config['Ra']
        Dec = config['Dec']
        rtot = []
        for season in range(1, config['nseasons']+1):
            mjd_min = config['MJD_min'] + float(season-1)*inter_season_gap
            mjd_max = mjd_min+config['season_length']

            for i, band in enumerate(bands):
                mjd = np.arange(mjd_min, mjd_max+cadence[band], cadence[band])
                mjd += shift_days[band]
                m5_coadded = self.m5_coadd(m5[band],
                                           Nvisits[band],
                                           Exposure_Time[band])
                myarr = np.array(mjd, dtype=[(self.mjdCol, 'f8')])
                myarr = rf.append_fields(myarr, [self.RaCol, self.DecCol, self.filterCol], [
                                         [Ra]*len(myarr), [Dec]*len(myarr), [band]*len(myarr)])
                myarr = rf.append_fields(myarr, [self.m5Col, self.nexpCol, self.exptimeCol, self.seasonCol], [
                                         [m5_coadded]*len(myarr), [Nvisits[band]]*len(myarr), [Nvisits[band]*Exposure_Time[band]]*len(myarr), [season]*len(myarr)])
                rtot.append(myarr)

        res = np.copy(np.concatenate(rtot))
        res.sort(order=self.mjdCol)

        self.Observations = res

    def m5_coadd(self, m5, Nvisits, Tvisit):
        m5_coadd = m5+1.25*np.log10(float(Nvisits)*Tvisit/30.)
        return m5_coadd


class TemplateData(object):
    """
    class to load template LC
    """

    def __init__(self, filename, band):
        """
        """
        self.fi = filename
        self.refdata = self.Stack()
        self.telescope = Telescope(airmass=1.1)
        self.blue_cutoff = 300.
        self.red_cutoff = 800.
        self.param_Fisher = ['X0', 'X1', 'Color']
        self.method = 'cubic'

        self.band = band
        idx = self.refdata['band'] == band
        lc_ref = self.refdata[idx]
        # load reference values (only once)
        phase_ref = lc_ref['phase']
        z_ref = lc_ref['z']
        flux_ref = lc_ref['flux']
        fluxerr_ref = lc_ref['fluxerr']
        self.gamma_ref = lc_ref['gamma'][0]
        self.m5_ref = np.unique(lc_ref['m5'])[0]
        self.dflux_interp = {}
        self.flux_interp = CloughTocher2DInterpolator(
            (phase_ref, z_ref), flux_ref, fill_value=1.e-5)
        self.fluxerr_interp = CloughTocher2DInterpolator(
            (phase_ref, z_ref), fluxerr_ref, fill_value=1.e-8)
        for val in self.param_Fisher:
            dflux = lc_ref['d'+val]
            self.dflux_interp[val] = CloughTocher2DInterpolator(
                (phase_ref, z_ref), dflux, fill_value=1.e-8)
            """
            dFlux[val] = griddata((phase, z), dflux, (phase_obs, yi_arr),
                                  method=self.method, fill_value=0.)
            """
        # this is to convert mag to flux in e per sec
        self.mag_to_flux_e_sec = {}
        mag_range = np.arange(14., 32., 0.1)
        for band in 'grizy':
            fluxes_e_sec = self.telescope.mag_to_flux_e_sec(
                mag_range, [band]*len(mag_range), [30]*len(mag_range))
            self.mag_to_flux_e_sec[band] = interpolate.interp1d(
                mag_range, fluxes_e_sec[:, 1], fill_value=0., bounds_error=False)

    def Stack(self):

        tab_tot = None
        f = h5py.File(self.fi, 'r')
        keys = f.keys()

        for kk in keys:

            tab_b = Table.read(self.fi, path=kk)
            if tab_tot is None:
                tab_tot = tab_b
            else:
                tab_tot = vstack([tab_tot, tab_b], metadata_conflicts='silent')

        return tab_tot

    def Fluxes(self, mjd_obs, param):

        z = param['z']
        daymax = param['DayMax']

        # observations (mjd, daymax, z) where to get the fluxes
        phase_obs = mjd_obs-daymax[:, np.newaxis]
        phase_obs = phase_obs/(1.+z[:, np.newaxis])  # phases of LC points
        z_arr = np.ones_like(phase_obs)*z[:, np.newaxis]

        flux = self.flux_interp((phase_obs, z_arr))

        return flux

    def Simulation(self, mjd_obs, m5_obs, exptime_obs, param):

        z = param['z']
        daymax = param['DayMax']

        # observations (mjd, daymax, z) where to get the fluxes
        phase_obs = mjd_obs-daymax[:, np.newaxis]
        phase_obs = phase_obs/(1.+z[:, np.newaxis])  # phases of LC points
        z_arr = np.ones_like(phase_obs)*z[:, np.newaxis]

        flux = self.flux_interp((phase_obs, z_arr))
        fluxerr = self.fluxerr_interp((phase_obs, z_arr))

        flux[flux <= 0] = 1.e-5
        """
        fluxerr_corr = self.FluxErrCorr(
            flux, m5_obs, exptime_obs, self.gamma_ref, self.m5_ref)
        
        fluxerr /= fluxerr_corr
        """
        tab = self.SelectSave(param, flux, fluxerr,
                              phase_obs, mjd_obs)
        return tab

    def FluxErrCorr(self, fluxes_obs, m5_obs, exptime_obs, gamma_ref, m5_ref):

        # Correct fluxes_err (m5 in generation probably different from m5 obs)
        gamma_obs = self.telescope.gamma(
            m5_obs, [self.band]*len(m5_obs), exptime_obs)
        mag_obs = -2.5*np.log10(fluxes_obs/3631.)

        gamma_tile = np.tile(gamma_obs, (len(mag_obs), 1))
        m5_tile = np.tile(m5_obs, (len(mag_obs), 1))
        srand_obs = self.Srand(gamma_tile, mag_obs, m5_tile)

        # srand_obs = self.srand(gamma_obs, mag_obs, m5_obs)
        # print('yes', band, m5_ref, gamma_ref, gamma_obs, mag_obs, srand_obs)

        m5 = np.asarray([m5_ref]*len(m5_obs))
        gamma = np.asarray([gamma_ref]*len(m5_obs))
        srand_ref = self.Srand(
            np.tile(gamma, (len(mag_obs), 1)), mag_obs, np.tile(m5, (len(mag_obs), 1)))

        correct_m5 = srand_ref/srand_obs
        return correct_m5

    def Srand(self, gamma, mag, m5):
        x = 10**(0.4*(mag-m5))
        return np.sqrt((0.04-gamma)*x+gamma*x**2)

    def FisherValues(self, param, phase_obs, fluxerr_obs):
        """
        idx = self.refdata['band'] == band
        lc_ref = self.refdata[idx]
        phase = lc_ref['phase']
        z = lc_ref['z']
        """
        yi_arr = np.ones_like(phase_obs)*param['z'][:, np.newaxis]
        dFlux = {}
        for val in self.param_Fisher:
            """
            dflux = lc_ref['d'+val]
            dFlux[val] = griddata((phase, z), dflux, (phase_obs, yi_arr),
                                  method=self.method, fill_value=0.)
            """
            dFlux[val] = self.dflux_interp[val]((phase_obs, yi_arr))

        FisherEl = {}
        for ia, vala in enumerate(self.param_Fisher):
            for jb, valb in enumerate(self.param_Fisher):
                if jb >= ia:
                    FisherEl[vala+valb] = dFlux[vala] * \
                        dFlux[valb]/fluxerr_obs**2
        return FisherEl

    def SelectSave(self, param, flux, fluxerr, phase, mjd):

        min_rf_phase = param['min_rf_phase']
        max_rf_phase = param['max_rf_phase']
        z = param['z']
        daymax = param['DayMax']

        # estimate element for Fisher matrix
        #FisherEl = self.FisherValues(param, phase, fluxerr)

        # flag for LC points outside the restframe phase range
        min_rf_phase = min_rf_phase[:, np.newaxis]
        max_rf_phase = max_rf_phase[:, np.newaxis]
        flag = (phase >= min_rf_phase) & (phase <= max_rf_phase)

        # flag for LC points outside the (blue-red) range
        mean_restframe_wavelength = np.array(
            [self.telescope.mean_wavelength[self.band]]*len(mjd))
        mean_restframe_wavelength = np.tile(
            mean_restframe_wavelength, (len(z), 1))/(1.+z[:, np.newaxis])
        flag &= (mean_restframe_wavelength > self.blue_cutoff) & (
            mean_restframe_wavelength < self.red_cutoff)
        flag_idx = np.argwhere(flag)

        # Now apply the flags to grab only interested values
        fluxes = np.ma.array(flux, mask=~flag)
        fluxes_err = np.ma.array(fluxerr, mask=~flag)
        mag = -2.5*np.log10(fluxes/3631.)
        phases = np.ma.array(phase, mask=~flag)
        snr_m5 = np.ma.array(flux/fluxerr, mask=~flag)
        obs_time = np.ma.array(
            np.tile(mjd, (len(mag), 1)), mask=~flag)
        """
        seasons = np.ma.array(
            np.tile(season, (len(mag_obs), 1)), mask=~flag)
        """
        z_vals = z[flag_idx[:, 0]]
        DayMax_vals = daymax[flag_idx[:, 0]]

        # Results are stored in an astropy Table
        tab = Table()
        tab.add_column(Column(fluxes[~fluxes.mask], name='flux'))
        tab.add_column(Column(fluxes_err[~fluxes_err.mask], name='fluxerr'))
        tab.add_column(Column(phases[~phases.mask], name='phase'))
        tab.add_column(Column(snr_m5[~snr_m5.mask], name='snr_m5'))
        tab.add_column(Column(mag[~mag.mask], name='mag'))
        tab.add_column(
            Column((2.5/np.log(10.))/snr_m5[~snr_m5.mask], name='magerr'))
        tab.add_column(Column(obs_time[~obs_time.mask], name='time'))

        tab.add_column(
            Column(['LSST::'+self.band]*len(tab), name='band',
                   dtype=h5py.special_dtype(vlen=str)))

        tab.add_column(Column([2.5*np.log10(3631)]*len(tab),
                              name='zp'))

        tab.add_column(
            Column(['ab']*len(tab), name='zpsys',
                   dtype=h5py.special_dtype(vlen=str)))

        # tab.add_column(Column(seasons[~seasons.mask], name='season'))
        tab.add_column(Column(z_vals, name='z'))
        tab.add_column(Column(DayMax_vals, name='DayMax'))
        """
        for key, vals in FisherEl.items():
            matel = np.ma.array(vals, mask=~flag)
            tab.add_column(Column(matel[~matel.mask], name='F_'+key))
        """
        return tab
