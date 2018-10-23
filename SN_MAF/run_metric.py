import matplotlib.pyplot as plt
import lsst.sims.maf.metricBundles as metricBundles
#import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.db as db
import lsst.sims.maf.utils as utils
from SN_Simulation import SNSimulation

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
