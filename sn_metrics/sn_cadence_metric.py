import numpy as np
from lsst.sims.maf.metrics import BaseMetric
from coadd_stacker import CoaddStacker
import healpy as hp
import numpy.lib.recfunctions as rf


class SNMetric(BaseMetric):
    """
    Measure mean m5, cadence and zlimit for SN
    """

    def __init__(self, metricName='SNMetric',
                 mjdCol='observationStartMJD', RaCol='fieldRA', DecCol='fieldDec',
                 filterCol='filter', m5Col='fiveSigmaDepth', exptimeCol='visitExposureTime',
                 nightCol='night', obsidCol='observationId', nexpCol='numExposures',
                 vistimeCol='visitTime', coadd=True, lim_sn=None, names_ref=None,
                 uniqueBlocks=False, config=None, **kwargs):

        self.mjdCol = mjdCol
        self.m5Col = m5Col
        self.filterCol = filterCol
        self.RaCol = RaCol
        self.DecCol = DecCol
        self.exptimeCol = exptimeCol
        self.seasonCol = 'season'
        self.nightCol = nightCol
        self.obsidCol = obsidCol
        self.nexpCol = nexpCol
        self.vistimeCol = vistimeCol

        cols = [self.nightCol, self.m5Col, self.filterCol, self.mjdCol, self.obsidCol,
                self.nexpCol, self.vistimeCol, self.exptimeCol, self.seasonCol]
        if coadd:
            cols += ['coadd']
        super(SNMetric, self).__init__(
            col=cols, metricDtype='object', metricName=metricName, **kwargs)

        self.filterNames = np.array(['u', 'g', 'r', 'i', 'z', 'y'])
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
        area = 9.6  # survey_area in sqdeg - 9.6 by default for DD
        if self.field_type == 'WFD':
            # in that case the survey area is the healpix area
            area = hp.nside2pixarea(
                config['Pixelisation']['nside'], degrees=True)

        # Load the reference Li file

        #self.Li = np.load(config['Reference File'])
        self.lim_sn = lim_sn
        self.names_ref = names_ref
        self.mag_range = config['Observations']['mag_range']
        self.dt_range = config['Observations']['dt_range']

    def run(self, dataSlice, slicePoint=None):
        # Cut down to only include filters in correct wave range.

        # print("attention",type(dataSlice),dataSlice.dtype)
        goodFilters = np.in1d(dataSlice['filter'], self.filterNames)
        dataSlice = dataSlice[goodFilters]
        if dataSlice.size == 0:
            return None
        dataSlice.sort(order=self.mjdCol)
        # print('dataslice',np.unique(dataSlice[['fieldRA','fieldDec','season','filter']]),dataSlice.dtype)
        time = dataSlice[self.mjdCol]-dataSlice[self.mjdCol].min()
        r = []
        fieldRA = np.mean(dataSlice[self.RaCol])
        fieldDec = np.mean(dataSlice[self.DecCol])
        band = np.unique(dataSlice[self.filterCol])[0]
        # for season  in np.unique(dataSlice[self.seasonCol]):
        seasons = [self.season]
        if self.season == -1:
            seasons = np.unique(dataSlice[self.seasonCol])
        for season in seasons:
            idx = dataSlice[self.seasonCol] == season
            sel = dataSlice[idx]
            bins = np.arange(np.floor(sel[self.mjdCol].min()), np.ceil(
                sel[self.mjdCol].max()), 1.)
            c, b = np.histogram(sel[self.mjdCol], bins=bins)
            cadence = 1. / c.mean()
            #time_diff = sel[self.mjdCol][1:]-sel[self.mjdCol][:-1]
            r.append((fieldRA, fieldDec, season, band,
                      np.mean(sel[self.m5Col]), cadence))
        # print(self.Li)

        res = np.rec.fromrecords(
            r, names=['fieldRA', 'fieldDec', 'season', 'band', 'm5_mean', 'cadence_mean'])

        idx = (res['m5_mean'] >= self.mag_range[0]) & (
            res['m5_mean'] <= self.mag_range[1])
        idx &= (res['cadence_mean'] >= self.dt_range[0]) & (
            res['cadence_mean'] <= self.dt_range[1])
        res = res[idx]
        # print(len(res))
        if len(res) > 0:
            if self.lim_sn is not None:
                for io, interp in enumerate(self.names_ref):
                    #zlims = [interp(xi, yi)[0] for xi, yi in zip(res['m5_mean'],res['cadence_mean'])]
                    zlims = self.lim_sn.Interp_griddata(io, res)
                    # print('interp',type(zlims))
                    zlims[np.isnan(zlims)] = -1
                    # print(io,zlims)
                    res = rf.append_fields(
                        res, 'zlim_'+self.names_ref[io], zlims)

        return res
