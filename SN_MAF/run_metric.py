import matplotlib.pyplot as plt
import lsst.sims.maf.metricBundles as metricBundles
#import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.db as db
import lsst.sims.maf.utils as utils
from SN_Simulation import SNSimulation
import argparse
import time
import yaml

parser = argparse.ArgumentParser(
    description='Run a SN metric from a configuration file')
parser.add_argument('config_filename',
                    help='Configuration file in YAML format.')

def run(config_filename):
    # YAML input file.
    config = yaml.load(open(config_filename))
    print(config)
    outDir ='Test'

    dbFile = '/data/pgris/sims_operation/Run_OpSim/enigma_1189_sqlite.db'
    dbFile='../flatiron/maf_local/sims_maf_contrib/tutorials/baseline2018a.db'
    #dbFile = '/data/pgris/sims_operation/Run_OpSim/clrlsstsrv_1068_sqlite.db'
    #opsimdb = utils.connectOpsimDb(dbFile)
    #resultsDb = db.ResultsDb(outDir=outDir)
    opsimdb = db.OpsimDatabase(dbFile)
    version = opsimdb.opsimVersion
    propinfo, proptags = opsimdb.fetchPropInfo()
    print('hello',proptags,propinfo)
    
    field='DD'
    #field='WFD'
    
    metric=SNSimulation(m5Col='fiveSigmaDepth', redshift=0.1, resolution=80.,Nbetween=10,singleDepthLimit=21.,peakGap=200.)
    #slicer = slicers.HealpixSlicer(nside=256)
    slicer=slicers.OpsimFieldSlicer()
    sqlconstraint = opsimdb.createSQLWhere(field, proptags)
    #sqlconstraint = utils.createSQLWhere('WFD', proptags)
    mb = metricBundles.MetricBundle(metric, slicer, sqlconstraint)

    mbD = {0:mb}

    resultsDb = db.ResultsDb(outDir=outDir)
    mbg =  metricBundles.MetricBundleGroup(mbD, opsimdb,
                                      outDir=outDir, resultsDb=resultsDb)

    mbg.runAll()

    #mbg.plotAll(closefigs=False)
    #plt.show()

def main(args):
    print('running')
    time_ref = time.time()
    run(args.config_filename)
    print('Time', time.time()-time_ref)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
