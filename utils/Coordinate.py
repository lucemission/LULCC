# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from index import YEAR_LATEST,EQ_BIOMASS,EQ_SOILSLOW,EX_BIOMASS,EX_SOILSLOW,EX_SOILRAPID,EX_PRODUCT1,EX_PRODUCT10,EX_PRODUCT100,EX_ATMOSPHERE, VIRGIN, SECOND, PASTURE, CROP, CLEARING, HARVEST, ABANDON, OTHERS, CLEAR_P1, CLEAR_P10, CLEAR_P100, CLEAR_SOIL, CLEAR_SRV, CLEAR_SRS, FORESTs, HARVEST_P1, HARVEST_P10, HARVEST_P100, HARVEST_SOIL_VIRGIN, HARVEST_SOIL_SECOND, HARVEST_MIN_SOIL
import Constants as c


class Coordinate:
    POOL_NUM = 7
    COVER_NUM = 4
    TRANS_NUM = 4
    PFT_NUM = 11
    FOREST_IDX = 6
    carbonEqu = np.zeros((1,2,COVER_NUM,PFT_NUM))
    # 
    carbonExcess = np.ones((1,POOL_NUM,COVER_NUM,TRANS_NUM,PFT_NUM))
    # time
    # Index: 0  1   2   3   4    5     6
    # Pool : B  SS  SR  P1  P10  P100  A
    # Cover: v  s   p   c
    
    # area = 30_802_500
    cover_pft_fraction  = np.zeros((4,PFT_NUM))
    cover_fraction      = np.zeros(4)
    # transition_pft_fraction = np.zeros((4,4,11))
    
    def __init__(self, start_year, transition, harvest_transition, initial_states, pft, duration, area):
        self.start_year = start_year
        self.transition_fraction = transition
        self.harvest_transition = harvest_transition
        self.states = initial_states
        self.pft = pft
        self.duration = duration
        self.area = area
    
    def visualize(self):
        
        def carbonOverTheYears(coorDict, coor, pool, cover, pft):
            if pool < 2:
                return coorDict[coor].carbonEqu[:, pool, cover, pft] + coorDict[coor].carbonExcess[:, pool, cover, :, pft].sum(axis=1)
            elif pool == 3:
                return coorDict[coor].carbonExcess[:, pool:pool+3, cover, :, pft].sum(axis=1).sum(axis=1)
            else:
                return coorDict[coor].carbonExcess[:, pool, cover, :, pft].sum(axis=1)
        carbonPools = dict()
        pools = ['biomass', 'slow soil', 'rapid soil', 'product', 'product 10', 'product 100', 'atm']
        colors = ['green', 'red', 'pink', 'tab:olive', 'tab:olive','tab:olive', 'blue']
        
        fig, ax = plt.subplots()
        
        POOL_NUM = 7
        cover = 1
        pft = 3
        
        years = np.arange(0, 101)
        for pool_idx in range(POOL_NUM):
            #     carbonPools[pools[pool_idx]] = 
            if (pool_idx == 4) or (pool_idx == 5):
                pass
            elif (pool_idx < 2):
                ax.plot(years, self.carbonEqu[:, pool_idx, cover, pft] + self.carbonExcess[:, pool_idx, cover, :, pft].sum(axis=1), label=pools[pool_idx], c=colors[pool_idx])
            elif (pool_idx == 3):
                ax.plot(years, self.carbonExcess[:, pool_idx:pool_idx+3, cover, :, pft].sum(axis=1).sum(axis=1), label=pools[pool_idx], c=colors[pool_idx])
            else: 
                ax.plot(years, self.carbonExcess[:, pool_idx, cover, :, pft].sum(axis=1), label=pools[pool_idx], c=colors[pool_idx])
        plt.axvline(x=10, color='grey', linestyle='--')
        plt.axvline(x=25, color='grey', linestyle='--')
        plt.axvline(x=40, color='grey', linestyle='--')
        plt.legend()
        plt.show()
        
    
    def run(self):
        for year in range(self.duration+1):
            self.updateFractions(year)
            self.event(year)
            self.relaxation()
        self.visualize()
        
    def updateFractions(self, year):
        if year == 0:
            self.carbonEqu[year,EQ_BIOMASS] = self.area*np.outer(self.states[year], self.pft[year])*c.carbonEquDense[EQ_BIOMASS]
        else:
            self.carbonEqu = np.append(self.carbonEqu, self.carbonEqu[YEAR_LATEST].copy()[np.newaxis, :], axis=0)
            self.carbonExcess = np.append(self.carbonExcess, self.carbonExcess[YEAR_LATEST].copy()[np.newaxis, :], axis=0)
        
        temp_carbonEquDense = c.carbonEquDense[EQ_BIOMASS].copy()
        temp_carbonEquDense[temp_carbonEquDense==0] = 1 # PLACEHOLDER
        
        # RUMUS 7
        self.cover_pft_fraction = self.carbonEqu[YEAR_LATEST,EQ_BIOMASS] / (self.area * temp_carbonEquDense)
        
        # RUMUS 8
        self.cover_fraction = self.cover_pft_fraction.sum(axis=1)
        
        # RUMUS 9
        self.temp_cover_fraction = self.cover_fraction.copy()
        self.temp_cover_fraction[self.temp_cover_fraction==0] = 1  
        self.transition_pft_fraction = np.tile(self.cover_pft_fraction[:,np.newaxis,:], (1,self.COVER_NUM,1)) * np.tile(self.transition_fraction[year][:,:,np.newaxis], (1,1,self.PFT_NUM)) / np.tile(self.temp_cover_fraction[:,np.newaxis,np.newaxis], (1,self.COVER_NUM,self.PFT_NUM))
        
    def event(self, year):
        
        # RUMUS 10
        self.BetaEquPrev = self.area * self.transition_pft_fraction * np.tile(c.carbonEquDense[EQ_BIOMASS][:,np.newaxis,:], (1,self.COVER_NUM,1))
        self.carbonEqu[YEAR_LATEST,EQ_BIOMASS] = self.carbonEqu[YEAR_LATEST,EQ_BIOMASS] - self.BetaEquPrev.sum(axis=1)
        
        # RUMUS 11 
        # (TIDAK DIPAKAI OLEH HARVEST)
        temp_cover_pft_fraction = self.cover_pft_fraction[:,np.newaxis,:].copy()
        temp_cover_pft_fraction[temp_cover_pft_fraction==0] = 1
        
        # 11
        # self.carbonExcess11 = self.carbonExcess[-1,0,:,0,:] - self.BetaExcessPrev[:,0,:]
        
        self.BetaExcessPrev = (np.tile(self.carbonExcess[YEAR_LATEST,EX_BIOMASS,:,np.newaxis,:,:], (1,self.COVER_NUM,1,1)) * self.createAltK() / temp_cover_pft_fraction[:,np.newaxis,:]) #.sum(axis=1)
        # APPEND
        self.carbonExcess[YEAR_LATEST,EX_BIOMASS,:,CLEARING,:] = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,:,CLEARING,:] -  self.BetaExcessPrev.sum(axis=1)[:,CLEARING,:]
        self.carbonExcess[YEAR_LATEST,EX_BIOMASS,:,ABANDON:,:] = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,:,ABANDON:,:] -  self.BetaExcessPrev.sum(axis=1)[:,ABANDON:,:]
        # self.BetaExcessPrev = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,:,:,:] * np.tile(self.transition_pft_fraction.sum(axis=1)[:,np.newaxis,:], (1,self.TRANS_NUM,1)) / temp_cover_pft_fraction #.sum(axis=1)
        # self.carbonExcess[YEAR_LATEST,EX_BIOMASS,:,:,:] = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,:,:,:] - self.BetaExcessPrev
        
        
        # RUMUS 12
        self.BetaExcessPrevSum = self.BetaExcessPrev.sum(axis=2)
        
        # # RUMUS 13
        self.Beta = self.BetaEquPrev + self.BetaExcessPrevSum
        
        # # RUMUS 14
        self.BetaEquNext = self.area * self.transition_pft_fraction.sum(axis=0) * c.carbonEquDense[EQ_BIOMASS]
        # APPEND
        self.carbonEqu[YEAR_LATEST,EQ_BIOMASS] = self.carbonEqu[YEAR_LATEST,EQ_BIOMASS] + self.BetaEquNext
        
        # # RUMUS 15
        # PERLU DIPERHATIKAN LAGI
        self.carbonExcess[YEAR_LATEST,EX_BIOMASS,:,CLEARING,:] = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,:,CLEARING,:] - self.BetaEquNext
        
        # # RUMUS 16
        self.carbonExcess[YEAR_LATEST,EX_PRODUCT1:EX_PRODUCT100+1,PASTURE:,CLEARING,:] = self.carbonExcess[YEAR_LATEST,EX_PRODUCT1:EX_PRODUCT100+1,PASTURE:,CLEARING,:] + np.tile(self.Beta[np.newaxis,:2,PASTURE:,:].sum(axis=1), (3,1,1)) * np.tile(c.clearing[:CLEAR_P100+1,np.newaxis,:], (1,2,1))
        self.carbonExcess[YEAR_LATEST,EX_SOILRAPID,PASTURE:,CLEARING,:] = self.carbonExcess[YEAR_LATEST,EX_SOILRAPID,PASTURE:,CLEARING,:] + (self.Beta[:SECOND+1,PASTURE:,:] * np.tile(c.clearing[CLEAR_SRV:,np.newaxis,:], (1,2,1))).sum(axis=0) * np.tile(c.clearing[np.newaxis,CLEAR_SOIL], (2,1))
        self.carbonExcess[YEAR_LATEST,EX_SOILSLOW,PASTURE:,CLEARING,:]  = self.carbonExcess[YEAR_LATEST,EX_SOILSLOW,PASTURE:,CLEARING,:]  + (self.Beta[:SECOND+1,PASTURE:,:] * np.tile(c.clearing[CLEAR_SRV:,np.newaxis,:], (1,2,1))).sum(axis=0) * (1 - np.tile(c.clearing[np.newaxis,CLEAR_SOIL], (2,1)))
        
        # # RUMUS 17
        self.SigmaEquPrev = self.area * self.transition_pft_fraction * np.tile(c.carbonEquDense[EQ_SOILSLOW][:,np.newaxis,:], (1,4,1))
        # APPEND
        self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW] = self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW] - self.SigmaEquPrev.sum(axis=1)
        
        
        # RUMUS 18
        self.SigmaExcessPrev = (np.tile(self.carbonExcess[YEAR_LATEST,EX_SOILSLOW,:,np.newaxis,:,:], (1,self.COVER_NUM,1,1)) * self.createAltK() / temp_cover_pft_fraction[:,np.newaxis,:]) #.sum(axis=1)
        self.carbonExcess[YEAR_LATEST,EX_SOILSLOW,:,CLEARING,:] = self.carbonExcess[YEAR_LATEST,EX_SOILSLOW,:,CLEARING,:] -  self.SigmaExcessPrev.sum(axis=1)[:,CLEARING,:]
        self.carbonExcess[YEAR_LATEST,EX_SOILSLOW,:,ABANDON:,:] = self.carbonExcess[YEAR_LATEST,EX_SOILSLOW,:,ABANDON:,:] -  self.SigmaExcessPrev.sum(axis=1)[:,ABANDON:,:]
        
        # RUMUS 19 
        self.SigmaExcessPrevSum = self.SigmaExcessPrev.sum(axis=2)
        
        # RUMUS 20
        self.Sigma = self.SigmaEquPrev + self.SigmaExcessPrevSum
        
        # # RUMUS 21
        self.SigmaEquNext = self.area * self.transition_pft_fraction * np.tile(c.carbonEquDense[EQ_SOILSLOW][np.newaxis,:,:], (self.COVER_NUM,1,1))
        self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW] = self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW] + self.SigmaEquNext.sum(axis=0)
        
        # RUMUS 22
        self.carbonExcess[YEAR_LATEST,EX_SOILRAPID,PASTURE:,CLEARING,:] = self.carbonExcess[YEAR_LATEST,EX_SOILRAPID,PASTURE:,CLEARING,:] + ((self.Sigma[:SECOND+1,PASTURE:,:] - self.SigmaEquNext[:SECOND+1,PASTURE:,:]) * np.tile(c.clearing[CLEAR_SRV:,np.newaxis,:], (1,2,1))).sum(axis=0)
        self.carbonExcess[YEAR_LATEST,EX_SOILSLOW,PASTURE:,CLEARING,:]  = self.carbonExcess[YEAR_LATEST,EX_SOILSLOW,PASTURE:,CLEARING,:]  + ((self.Sigma[:SECOND+1,PASTURE:,:] - self.SigmaEquNext[:SECOND+1,PASTURE:,:]) * (1 - np.tile(c.clearing[CLEAR_SRV:,np.newaxis,:], (1,2,1)))).sum(axis=0)
        
        # RUMUS 23
        # Abandonment
        self.carbonExcess[YEAR_LATEST,EX_BIOMASS,SECOND,ABANDON,:] = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,SECOND,ABANDON,:]  + (self.Beta[PASTURE:,SECOND,:].sum(axis=0) - self.BetaEquNext[SECOND,:])
        # Others
        self.carbonExcess[YEAR_LATEST,EX_BIOMASS,PASTURE,OTHERS,:] = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,PASTURE,OTHERS,:]  + (self.Beta[CROP,PASTURE,:].sum(axis=0) - self.BetaEquNext[PASTURE,:])
        self.carbonExcess[YEAR_LATEST,EX_BIOMASS,CROP,OTHERS,:] = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,CROP,OTHERS,:]  + (self.Beta[PASTURE,CROP,:].sum(axis=0) - self.BetaEquNext[CROP,:])
        
        # RUMUS 24
        # Abandonment
        self.carbonExcess[YEAR_LATEST,EX_BIOMASS,SECOND,ABANDON,:]  = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,SECOND,ABANDON,:]  + (self.SigmaEquPrev[PASTURE:,SECOND,:] - self.SigmaEquNext[PASTURE:,SECOND,:]).sum(axis=0)
        # Others
        self.carbonExcess[YEAR_LATEST,EX_BIOMASS,PASTURE,OTHERS,:]  = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,PASTURE,OTHERS,:]  + (self.SigmaEquPrev[CROP,PASTURE,:] - self.SigmaEquNext[CROP,PASTURE,:]).sum(axis=0)
        self.carbonExcess[YEAR_LATEST,EX_BIOMASS,CROP,OTHERS,:]  = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,CROP,OTHERS,:]  + (self.SigmaEquPrev[PASTURE,CROP,:] - self.SigmaEquNext[PASTURE,CROP,:]).sum(axis=0)
        
        
        # # Prelude rumus 25
        
        temp_harvest_cover_fraction = self.createHarvestCoverDivisor()
        
        # Rumus 9 Harvest
        self.harvest_transition_pft_fraction = self.cover_pft_fraction[VIRGIN:SECOND+1,:] * self.harvest_transition[year] / temp_harvest_cover_fraction
        # # # Rumus 10 Harvest
        self.BetaHarvestEquPrev = self.area * self.harvest_transition_pft_fraction * c.carbonEquDense[EQ_BIOMASS,VIRGIN:SECOND+1]
        self.carbonEqu[YEAR_LATEST,EQ_BIOMASS,VIRGIN:SECOND+1,:] = self.carbonEqu[YEAR_LATEST,EQ_BIOMASS,VIRGIN:SECOND+1,:] - self.BetaHarvestEquPrev
        
        # # # RUMUS 25
        self.BetaHarvestExcessPrev = self.carbonExcess[YEAR_LATEST, EX_BIOMASS,VIRGIN:SECOND+1, HARVEST,:] * self.harvest_transition_pft_fraction / temp_harvest_cover_fraction
        
       # Rumus 14 Harvest
        self.BetaHarvestEquNext = self.area * self.harvest_transition_pft_fraction.sum(axis=0) * c.carbonEquDense[EQ_BIOMASS, SECOND]
        self.carbonEqu[YEAR_LATEST,EQ_BIOMASS, SECOND, :] = self.carbonEqu[YEAR_LATEST,EQ_BIOMASS, SECOND, :] + self.BetaHarvestEquNext
        
        # RUMUS 26
        self.carbonExcess[YEAR_LATEST, EX_BIOMASS, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, EX_BIOMASS,SECOND,HARVEST,:] - (self.BetaHarvestEquNext + self.BetaHarvestExcessPrev.sum(axis=0))
        
        # RUMUS 27
        self.BetaHarvest = self.BetaHarvestEquPrev + self.BetaHarvestExcessPrev
        
        # RUMUS 28
        self.carbonExcess[YEAR_LATEST, EX_SOILRAPID, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, EX_SOILRAPID, SECOND, HARVEST, :] + (self.BetaHarvest * c.harvest[HARVEST_SOIL_VIRGIN: HARVEST_SOIL_SECOND+1]).sum(axis=0)
        
        # RUMUS 29
        self.carbonExcess[YEAR_LATEST, EX_PRODUCT1, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, EX_PRODUCT1, SECOND, HARVEST, :] + (self.BetaHarvest * (1-c.harvest[HARVEST_SOIL_VIRGIN: HARVEST_SOIL_SECOND+1])).sum(axis=0) * c.harvest[HARVEST_P1]
        self.carbonExcess[YEAR_LATEST, EX_PRODUCT10, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, EX_PRODUCT10, SECOND, HARVEST, :] + (self.BetaHarvest * (1-c.harvest[HARVEST_SOIL_VIRGIN: HARVEST_SOIL_SECOND+1])).sum(axis=0) * c.harvest[HARVEST_P10]
        self.carbonExcess[YEAR_LATEST, EX_PRODUCT100, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, EX_PRODUCT100, SECOND, HARVEST, :] + (self.BetaHarvest * (1-c.harvest[HARVEST_SOIL_VIRGIN: HARVEST_SOIL_SECOND+1])).sum(axis=0) * c.harvest[HARVEST_P100]
        
        # Rumus 17 Harvest
        self.SigmaHarvestEquPrev = self.area * self.harvest_transition_pft_fraction * c.carbonEquDense[EQ_SOILSLOW,VIRGIN:SECOND+1]
        self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW,VIRGIN:SECOND+1,:] = self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW,VIRGIN:SECOND+1,:] - self.SigmaHarvestEquPrev
        
        # RUMUS 30
        self.SigmaHarvestExcessPrev = self.carbonExcess[YEAR_LATEST, EX_SOILSLOW,VIRGIN:SECOND+1, HARVEST,:] * self.harvest_transition_pft_fraction / temp_harvest_cover_fraction
        
        # RUMUS 31
        self.SigmaHarvest = self.SigmaHarvestEquPrev + self.SigmaHarvestExcessPrev
        
        # rumus 21 Harvest
        self.SigmaHarvestEquNext = self.area * self.harvest_transition_pft_fraction.sum(axis=0) * c.carbonEquDense[EQ_SOILSLOW, SECOND]
        self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW, SECOND, :] = self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW, SECOND, :] + self.SigmaHarvestEquNext
        
        # RUMUS 32
        # PERLU DIPERHATIKAN LAGI
        self.SigmaRapid = np.maximum(0, self.SigmaHarvest.sum(axis=0) - self.area * self.harvest_transition_pft_fraction.sum(axis=0) * c.harvestMinimumSoil)
        self.carbonExcess[YEAR_LATEST, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, SECOND, HARVEST, :] + self.SigmaRapid
        
        # RUMUS 33
        self.carbonExcess[YEAR_LATEST, EX_SOILSLOW, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, EX_SOILSLOW, SECOND, HARVEST, :] - self.SigmaRapid + self.SigmaHarvestEquPrev.sum(axis=0) - self.SigmaHarvestEquNext
        
    # RUMUS 6
    def relaxation(self):
        relaxed_carbon = np.zeros((self.POOL_NUM-1,self.COVER_NUM,self.TRANS_NUM,self.PFT_NUM))
        
        relaxed_carbon[:EX_SOILRAPID+1,:,:,:] = self.carbonExcess[YEAR_LATEST,:EX_SOILRAPID+1,:,:,:]*np.exp(-1/c.Timescale)
        relaxed_carbon[EX_PRODUCT1,:,:,:] = self.carbonExcess[YEAR_LATEST,EX_PRODUCT1,:,:,:]*np.exp(-1/np.tile(0.534, (4,4,11)))
        relaxed_carbon[EX_PRODUCT10,:,:,:] = self.carbonExcess[YEAR_LATEST,EX_PRODUCT10,:,:,:]*np.exp(-1/np.tile(5.34, (4,4,11)))
        relaxed_carbon[EX_PRODUCT100,:,:,:] = self.carbonExcess[YEAR_LATEST,EX_PRODUCT100,:,:,:]*np.exp(-1/np.tile(53.4, (4,4,11)))
        
    
        
        self.carbonExcess[YEAR_LATEST,EX_ATMOSPHERE,:,:,:] = self.carbonExcess[YEAR_LATEST,EX_ATMOSPHERE,:,:,:] + relaxed_carbon.sum(axis=0)
        self.carbonExcess[YEAR_LATEST,:EX_ATMOSPHERE,:,:,:] = self.carbonExcess[YEAR_LATEST,:EX_ATMOSPHERE,:,:,:] - relaxed_carbon
    
    def createAltK(self):
        result = np.zeros((self.transition_fraction.shape[1], self.transition_fraction.shape[2], self.transition_fraction.shape[2], self.cover_pft_fraction.shape[1]))
        
        # virgin - secondary - harvest
        result[0,1,1,:] = self.transition_pft_fraction[0,1,:]
        # virgin - pasture - clearing
        result[0,2,0,:] = self.transition_pft_fraction[0,2,:]
        # virgin - crop - clearing
        result[0,3,0,:] = self.transition_pft_fraction[0,3,:]
        
        # secondary - pasture - clearing
        result[1,2,0,:] = self.transition_pft_fraction[1,2,:]
        # secondary - crop - clearing
        result[1,3,0,:] = self.transition_pft_fraction[1,3,:]
        
        # pasture - secondary - abandonment
        result[2,1,2,:] = self.transition_pft_fraction[2,1,:]
        # pasture - crop - others
        result[2,3,3,:] = self.transition_pft_fraction[2,3,:]
        
        # crop - secondary - abandonment
        result[3,1,2,:] = self.transition_pft_fraction[3,1,:]
        # crop - pasture - others
        result[3,2,3,:] = self.transition_pft_fraction[3,2,:]
        return result
    
    def createHarvestCoverDivisor(self):
        divisor = self.cover_fraction.copy()
        divisor[divisor==0] = 1 
        divisor = np.tile(divisor[VIRGIN:SECOND+1,np.newaxis],(1,11))
        return divisor
    
    # def createHarvestCoverPftFraction(self):\
        
    #     result = np.zeros((2,2,self.FOREST_IDX))
        
    #     result[VIRGIN, FORESTs,:] = self.cover_pft_fraction[VIRGIN, 0:self.FOREST_IDX]
    #     result[VIRGIN, NONFOREST,:] = self.cover_pft_fraction[VIRGIN, self.FOREST_IDX:]
    #     # result[VIRGIN, NONFOREST, :-1] = self.cover_pft_fraction[VIRGIN, self.FOREST_IDX:]
    #     result[SECOND, FOREST,:] = self.cover_pft_fraction[SECOND, 0:self.FOREST_IDX]
    #     result[SECOND, NONFOREST,:] = self.cover_pft_fraction[SECOND, self.FOREST_IDX:]
    #     # result[SECOND, NONFOREST,:-1] = self.cover_pft_fraction[SECOND, self.FOREST_IDX:]
    #     return result
        
    def createTransitionK(self):
        result = np.zeros((self.transition_fraction.shape[1], self.transition_fraction.shape[2], self.cover_pft_fraction.shape[1]))
        
        # # virgin - secondary - harvest
        # result[0,1,1,:] = self.transition_pft_fraction[0,1,:]
        # # virgin - pasture - clearing
        # result[0,2,0,:] = self.transition_pft_fraction[0,2,:]
        # # virgin - crop - clearing
        # result[0,3,0,:] = self.transition_pft_fraction[0,3,:]
        
        # # secondary - pasture - clearing
        # result[1,2,0,:] = self.transition_pft_fraction[1,2,:]
        # # secondary - crop - clearing
        # result[1,3,0,:] = self.transition_pft_fraction[1,3,:]
        
        # # pasture - secondary - abandonment
        # result[2,1,2,:] = self.transition_pft_fraction[2,1,:]
        # # pasture - crop - others
        # result[2,3,3,:] = self.transition_pft_fraction[2,3,:]
        
        # # crop - secondary - abandonment
        # result[3,1,2,:] = self.transition_pft_fraction[3,1,:]
        # # crop - pasture - others
        # result[3,2,3,:] = self.transition_pft_fraction[3,2,:]
        
        # virgin clearing
        result[0,0,:] = self.transition_pft_fraction[0,2:,:].sum(axis=0)
        # virgin harvest
        result[0,1,:] = self.transition_pft_fraction[0,1,:]
        
        # secondary clearing
        result[1,0] = self.transition_pft_fraction[1,2:,:].sum(axis=0)
        # result[1,1] untuk secondary harvest
        
        # pasture abandonment
        result[2,2] = self.transition_pft_fraction[2,1,:]
        # pasture others
        result[2,3] = self.transition_pft_fraction[2,3,:]
        
        # crop abandonment
        result[3,2] = self.transition_pft_fraction[3,1,:]
        # crop others
        result[3,3] = self.transition_pft_fraction[3,2,:]
        
        return result
        
        
    
        
    
        
