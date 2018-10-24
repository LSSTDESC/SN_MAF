from lsst.sims.maf.stackers import BaseStacker
import numpy as np
import numpy.lib.recfunctions as rf

__all__ = ['CoaddStacker']

class CoaddStacker(BaseStacker):

    colsAdded = ['coadd']
    def __init__(self,RaCol='fieldRA',DecCol='fieldDec',m5Col='fiveSigmaDepth', nightcol = 'night',filterCol='filter',nightCol='night', numExposuresCol='numExposures', visitTimeCol='visitTime', visitExposureTimeCol='visitExposureTime'):
        self.colsReq = [RaCol, DecCol,m5Col,filterCol,nightCol,numExposuresCol,visitTimeCol,visitExposureTimeCol]
        self.RaCol = RaCol
        self.DecCol = DecCol
        self.nightCol = nightCol
        self.filterCol = filterCol
        self.m5Col = m5Col
        self.numExposuresCol = numExposuresCol
        self.visitTimeCol = visitTimeCol
        self.visitExposureTimeCol = visitExposureTimeCol
        
        self.units = ['int']
        
    def _run(self, simData, cols_present = False):
        if cols_present:
            # Column already present in data; assume it is correct and does not need recalculating.
            return simData
        self.dtype = simData.dtype
        r= []
        for ra, dec, band in np.unique(simData[[self.RaCol, self.DecCol, self.filterCol]]):
            idx = np.abs(simData[self.RaCol]-ra) < 1.e-5
            idx &= np.abs(simData[self.DecCol]-dec) < 1.e-5
            idx &= simData[self.filterCol] == band
            
            sel = simData[idx]
            #print('there',len(sel))
            for night in np.unique(sel[self.nightCol]):
                idxb = sel[self.nightCol] == night
                #print(ra,dec,band, night,len(sel[idxb]),sel[idxb][self.numExposuresCol],sel[idxb][self.visitExposureTimeCol],sel[idxb][self.visitTimeCol])
                r.append(tuple(self.Fill(sel[idxb])))

        #print(r)
        myarray = np.array(r, dtype = self.dtype)
        return myarray

    def Fill(self, tab):

        r = []

        for colname in self.dtype.names:
            if colname not in [self.m5Col, self.numExposuresCol, self.visitTimeCol, self.visitExposureTimeCol,self.filterCol]:
                if colname =='coadd':
                    r.append(1)
                else:
                    r.append(np.median(tab[colname]))
            if colname == self.m5Col:
                r.append(self.m5_coadd(tab[self.m5Col]))
            if colname in [self.numExposuresCol, self.visitTimeCol, self.visitExposureTimeCol]:
                r.append(np.sum(tab[colname]))
            if colname == self.filterCol:
                r.append(np.unique(tab[self.filterCol])[0])
       
        return r
    
    def m5_coadd(self,m5):

        fluxes=10**(-0.4*m5)
        sigmas=fluxes/5.
        sigma_tot=1./np.sqrt(np.sum(1./sigmas**2))
        flux_tot=5.*sigma_tot
    
        return -2.5*np.log10(flux_tot)
