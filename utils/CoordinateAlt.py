# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from index import YEAR_LATEST,EQ_BIOMASS,EQ_SOILSLOW,EX_BIOMASS,EX_SOILSLOW,EX_SOILRAPID,EX_PRODUCT1,EX_PRODUCT10,EX_PRODUCT100,EX_ATMOSPHERE, VIRGIN, SECOND, PASTURE, CROP, CLEARING, HARVEST, ABANDON, OTHERS, CLEAR_P1, CLEAR_P10, CLEAR_P100, CLEAR_SOIL, CLEAR_SRV, CLEAR_SRS, FORESTs, HARVEST_P1, HARVEST_P10, HARVEST_P100, HARVEST_SOIL_VIRGIN, HARVEST_SOIL_SECOND, HARVEST_MIN_SOIL
import equations
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
        self.area = np.squeeze(area)
    
    # def visualize(self):
    #     years = np.linspace()
    #     fig,ax = plt.subplots()
        
    
    def run(self):
        for year in range(self.duration+1):
            self.updateFractions(year)
            self.event(year)
            self.relaxation()
        
    def updateFractions(self, year):
        if year == 0:
            # Jika tahun pertama (ke-0), lakukan inisialisasi nilai carbonEqu
            self.carbonEqu[year,EQ_BIOMASS] = self.area*np.outer(self.states[year], self.pft[year])*c.carbonEquDense[EQ_BIOMASS]
        else:
            # Jika tidak, copy carbonEqu tahun lalu dan append pada variabel carbonEqu
            self.carbonEqu = np.append(self.carbonEqu, self.carbonEqu[YEAR_LATEST].copy()[np.newaxis, :], axis=0)
            # copy juga carbonExcess tahun lalu dan append pada variabel carbonExcess
            self.carbonExcess = np.append(self.carbonExcess, self.carbonExcess[YEAR_LATEST].copy()[np.newaxis, :], axis=0)
        
        temp_carbonEquDense = c.carbonEquDense[EQ_BIOMASS].copy()
        temp_carbonEquDense[temp_carbonEquDense==0] = 1 # PLACEHOLDER
        
        # RUMUS 7
        # self.cover_pft_fraction = self.carbonEqu[YEAR_LATEST,EQ_BIOMASS] / (self.area * temp_carbonEquDense)
        self.cover_pft_fraction = equations.equation7(self.carbonEqu[YEAR_LATEST,EQ_BIOMASS], self.area, temp_carbonEquDense)
        
        # RUMUS 8
        self.cover_fraction = self.cover_pft_fraction.sum(axis=1)
        
        # RUMUS 9
        self.temp_cover_fraction = self.cover_fraction.copy()
        self.temp_cover_fraction[self.temp_cover_fraction==0] = 1  

        self.transition_pft_fraction = np.tile(self.cover_pft_fraction[:,np.newaxis,:], (1,self.COVER_NUM,1)) * np.tile(self.transition_fraction[year][:,:,np.newaxis], (1,1,self.PFT_NUM)) / np.tile(self.temp_cover_fraction[:,np.newaxis,np.newaxis], (1,self.COVER_NUM,self.PFT_NUM))

        # RUMUS 9 HARVEST
        self.temp_harvest_cover_fraction = self.createHarvestCoverDivisor()
        self.harvest_transition_pft_fraction = self.cover_pft_fraction[VIRGIN:SECOND+1,:] * self.harvest_transition[year] / self.temp_harvest_cover_fraction
        
    def event(self, year):
        self.clearing()
        self.abandonment()
        self.others_pc()
        self.others_cp()
        self.harvest()
        
    def clearing(self):
        pasture = 0
        crop = 1

        sources = [VIRGIN, SECOND]
        targets = [pasture, crop]
        TARGETS = [PASTURE, CROP]

        self.biomass_clearing(sources, targets, TARGETS)
        self.soil_clearing(sources, targets, TARGETS)

    def abandonment(self):
        pasture = 0
        crop = 1
        second = 0

        sources  = [pasture, crop]
        SOURCES = [PASTURE, CROP]
        TARGET = SECOND
        target = second

        self.biomass_abandonment(sources, SOURCES, target, TARGET)
        self.soil_abandonment(sources, SOURCES, target, TARGET)
    
    def others_pc(self):
        source = pasture = 0
        target = crop = 0

        SOURCE = PASTURE
        TARGET = CROP

        self.biomass_others(source, SOURCE, target, TARGET)
        self.soil_others(source, SOURCE, target, TARGET)

    def others_cp(self):
        source = crop = 0
        target = pasture = 0

        SOURCE = CROP
        TARGET = PASTURE

        self.biomass_others(source, SOURCE, target, TARGET)
        self.soil_others(source, SOURCE, target, TARGET)

    def harvest(self):

        target = second = 0
        SOURCES = [VIRGIN, SECOND]
        TARGET = SECOND

        self.biomass_harvest(SOURCES, target, TARGET)
        self.soil_harvest(SOURCES, target, TARGET)

    # RUMUS 6
    def relaxation(self):
        relaxed_carbon = np.zeros((self.POOL_NUM-1,self.COVER_NUM,self.TRANS_NUM,self.PFT_NUM))
        
        relaxed_carbon[:EX_SOILRAPID+1,:,:,:] = self.carbonExcess[YEAR_LATEST,:EX_SOILRAPID+1,:,:,:]*np.exp(-1/c.Timescale)
        relaxed_carbon[EX_PRODUCT1,:,:,:] = self.carbonExcess[YEAR_LATEST,EX_PRODUCT1,:,:,:]*np.exp(-1/np.tile(0.534, (4,4,11)))
        relaxed_carbon[EX_PRODUCT10,:,:,:] = self.carbonExcess[YEAR_LATEST,EX_PRODUCT10,:,:,:]*np.exp(-1/np.tile(5.34, (4,4,11)))
        relaxed_carbon[EX_PRODUCT100,:,:,:] = self.carbonExcess[YEAR_LATEST,EX_PRODUCT100,:,:,:]*np.exp(-1/np.tile(53.4, (4,4,11)))
        
        self.carbonExcess[YEAR_LATEST,EX_ATMOSPHERE,:,:,:] = self.carbonExcess[YEAR_LATEST,EX_ATMOSPHERE,:,:,:] + relaxed_carbon.sum(axis=0)
        self.carbonExcess[YEAR_LATEST,:EX_ATMOSPHERE,:,:,:] = self.carbonExcess[YEAR_LATEST,:EX_ATMOSPHERE,:,:,:] - relaxed_carbon
    
    def equation10and17(self, pool, source, target):
        """
        Menghitung jumlah karbon equilibrium yang dilepaskan

        Parameters:
            pool (int) : indeks tipe pool karbon (Biomass, Slow Soil, Rapid Soil, 1-year Product, 10-year Product, 100-year Product, Atmosphere)
            source (int) : indeks source cover type absolut (virgin, secondary, pasture, crop)
            target (int) : indeks target cover type absolut (virgin, secondary, pasture, crop)

        returns
            jumlah karbon ekuilibrium yang dilepaskan dari source
        """
        # print('eq 10')
        # print(self.transition_pft_fraction[source, target, :].shape)
        symbolEquPrev = self.area * self.transition_pft_fraction[source, target, :] * c.carbonEquDense[pool, source,:]
        self.carbonEqu[YEAR_LATEST,pool, source, :] = self.carbonEqu[YEAR_LATEST,pool, source, :] - symbolEquPrev

        return symbolEquPrev

    def equation11and18(self, pool, source, target):
        """
        Menghitung lepasnya jumlah karbon excess yang dilepaskan

        Parameters:
            pool (int) : indeks tipe pool karbon (Biomass, Slow Soil, Rapid Soil, 1-year Product, 10-year Product, 100-year Product, Atmosphere)
            source (int) : indeks source cover type absolut (virgin, secondary, pasture, crop)
            target (int) : indeks target cover type absolut (virgin, secondary, pasture, crop)

        returns
            jumlah karbon excess yang dilepaskan dari source
        """
        temp_cover_pft_fraction = self.cover_pft_fraction[source,np.newaxis,:].copy()
        temp_cover_pft_fraction[temp_cover_pft_fraction==0] = 1
        temp_cover_pft_fraction = np.tile(temp_cover_pft_fraction, (self.TRANS_NUM, 1))

        
        BetaExcessPrev = self.carbonExcess[YEAR_LATEST,pool,source, :,:] * np.tile(self.transition_pft_fraction[source, target,np.newaxis, :], (self.TRANS_NUM,1)) / temp_cover_pft_fraction
        self.carbonExcess[YEAR_LATEST,pool,source, :,:] = self.carbonExcess[YEAR_LATEST,pool,source, :,:] - BetaExcessPrev
        return BetaExcessPrev

    def equation12and19(self, symbolExcessPrev, source, target):
        return symbolExcessPrev[source, target, :, :].sum(axis=0)

    def equation13and20(self,symbolEqu,symbolExcessPrevSum, source, target):
        return symbolEqu[source, target, :] +symbolExcessPrevSum[source, target, :]

    def equation14(self, sources, target):
        """
        Menghitung jumlah karbon ekuilibrium yang diterima oleh target cover type

        Parameters:
            target (int) : indeks target cover type absolut
        """
        if type(sources) == list:
            transition_pft = self.transition_pft_fraction[sources[0]:sources[1]+1].sum(axis=0)[target,:]
        elif type(sources) == int:
            transition_pft = self.transition_pft_fraction[sources,target,:]
        BetaEquNext = np.squeeze(self.area * transition_pft * c.carbonEquDense[EQ_BIOMASS, target, :])
        self.carbonEqu[YEAR_LATEST,EQ_BIOMASS, target, :] = self.carbonEqu[YEAR_LATEST,EQ_BIOMASS, target, :] + BetaEquNext
        return BetaEquNext

    def equation15(self, betaEquNext, target1, trans):
        self.carbonExcess[YEAR_LATEST, EX_BIOMASS, target1, trans, :] = self.carbonExcess[YEAR_LATEST, EX_BIOMASS, target1, trans, :] - betaEquNext

    def equation16_product(self, pool, target, beta, clear_frac):
        self.carbonExcess[YEAR_LATEST, pool, target, CLEARING, :] = self.carbonExcess[YEAR_LATEST, pool, target, CLEARING, :] + beta * c.clearing[clear_frac, :]

    def equation16_soil(self, pool, target, beta):
        if pool == EX_SOILRAPID:
            self.carbonExcess[YEAR_LATEST, pool, target, CLEARING, :] = self.carbonExcess[YEAR_LATEST, pool, target, CLEARING, :] + (beta * c.clearing[CLEAR_SRV:CLEAR_SRS+1,:]).sum(axis=0) * c.clearing[CLEAR_SOIL]
        elif pool == EX_SOILSLOW:
            self.carbonExcess[YEAR_LATEST, pool, target, CLEARING, :] = self.carbonExcess[YEAR_LATEST, pool, target, CLEARING, :] + (beta * 1-c.clearing[CLEAR_SRV:CLEAR_SRS+1,:]).sum(axis=0) * c.clearing[CLEAR_SOIL]

    def equation21(self, pool, SOURCE, TARGET):
        symbolEquNext = self.area * self.transition_pft_fraction[SOURCE, TARGET, :] * c.carbonEquDense[pool, TARGET,:]
        self.carbonEqu[YEAR_LATEST,pool, TARGET, :] = self.carbonEqu[YEAR_LATEST,pool, TARGET, :] + np.squeeze(symbolEquNext)
        return symbolEquNext

    def equation22(self, pool, target1, target2, trans, sigma, sigmaEquNext):
        if pool == EX_SOILRAPID:
            self.carbonExcess[YEAR_LATEST, pool, target1, trans, :] = self.carbonExcess[YEAR_LATEST, pool, target1, trans, :] + ((sigma[:, target2,:] - sigmaEquNext[:,target2,:])*c.clearing[CLEAR_SRV:CLEAR_SRS+1,:]).sum(axis=0)
        elif pool == EX_SOILSLOW:
            self.carbonExcess[YEAR_LATEST, pool, target1, trans, :] = self.carbonExcess[YEAR_LATEST, pool, target1, trans, :] + ((sigma[:, target2,:] - sigmaEquNext[:,target2,:])*(1-c.clearing[CLEAR_SRV:CLEAR_SRS+1,:])).sum(axis=0)

    def equation23(self, TARGET, trans, betaEquNext, beta):
        """
        Menambahkan karbon excess ke target cover type

        Parameters:
            TARGET (int) : indeks target cover type absolut
            target (int) : indeks target cover type kontektual
            trans (int) : indeks transisi
        """
        self.carbonExcess[YEAR_LATEST, EX_BIOMASS, TARGET, trans, :] = self.carbonExcess[YEAR_LATEST, EX_BIOMASS, TARGET, trans, :] - betaEquNext + np.squeeze(beta.sum(axis=0))

    def equation24(self, TARGET, trans, sigmaEquPrev, sigmaEquNext):
        self.carbonExcess[YEAR_LATEST, EX_SOILSLOW, TARGET, trans, :] = self.carbonExcess[YEAR_LATEST, EX_SOILSLOW, TARGET, trans, :] + np.squeeze(sigmaEquPrev.sum(axis=0)) - sigmaEquNext

    def equation25a(self, SOURCE, TARGET):
        temp_cover_pft_fraction = self.cover_pft_fraction[SOURCE,np.newaxis,:].copy()
        temp_cover_pft_fraction[temp_cover_pft_fraction==0] = 1
        temp_cover_pft_fraction = np.tile(temp_cover_pft_fraction, (self.TRANS_NUM, 1))

        BetaExcessPrev = self.carbonExcess[YEAR_LATEST,EX_BIOMASS,SOURCE, :,:] * np.tile(self.transition_pft_fraction[SOURCE, TARGET,np.newaxis, :], (self.TRANS_NUM,1)) / temp_cover_pft_fraction
        return BetaExcessPrev

    def equation26(self, betaEquNext, betaExcess):
        self.carbonExcess[YEAR_LATEST, EX_BIOMASS, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST,EX_BIOMASS, SECOND, HARVEST, :] - (betaEquNext + np.squeeze(betaExcess.sum(axis=0)))

    def equation27(self, betaEquPrev, betaExcessPrevSum):
        return betaEquPrev + betaExcessPrevSum
        
    def equation28(self, betaHarvest):
        self.carbonExcess[YEAR_LATEST, EX_SOILRAPID, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, EX_SOILRAPID, SECOND, HARVEST, :] + (betaHarvest * c.harvest[HARVEST_SOIL_VIRGIN: HARVEST_SOIL_SECOND+1]).sum(axis=0)

    def equation29(self, betaHarvest):
        remainder = 1 - c.harvest[HARVEST_SOIL_VIRGIN: HARVEST_SOIL_SECOND+1]
        products = [EX_PRODUCT1, EX_PRODUCT10, EX_PRODUCT100]
        product_constants = [HARVEST_P1, HARVEST_P10, HARVEST_P100]

        for i in range(len(products)):
            self.carbonExcess[YEAR_LATEST, products[i], SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, products[i], SECOND, HARVEST, :] + (betaHarvest * remainder).sum(axis=0) * c.harvest[product_constants[i]]

    def equation30a(self, source):
        temp_cover_pft_fraction = self.cover_pft_fraction[source,np.newaxis,:].copy()
        temp_cover_pft_fraction[temp_cover_pft_fraction==0] = 1
        temp_cover_pft_fraction = np.tile(temp_cover_pft_fraction, (self.TRANS_NUM, 1))
        
        return self.carbonExcess[YEAR_LATEST, EX_SOILSLOW, source,:,:] * np.tile(self.transition_pft_fraction[source, SECOND,np.newaxis:], (self.TRANS_NUM, 1)) / temp_cover_pft_fraction

    def biomass_clearing(self, sources, targets, TARGETS):
        
        # print()
        # print('clearing')
        
        clearBetaEquPrev, clearBetaExcessPrev, clearBetaExcessPrevSum, clearBeta = self.initializeSymbols(len(sources), len(targets), self.TRANS_NUM, self.PFT_NUM)

        for source in sources:
            for j in range(len(targets)):
                target = targets[j]
                TARGET = TARGETS[j]

                # clearBetaEquPrev[source, targets[j], :] = self.equation10and17(EQ_BIOMASS, source, TARGETS[j])
                clearBetaEquPrev[source, target, :] = equations.equation10and17(self.area, self.transition_pft_fraction[source, target], c.carbonEquDense[EQ_BIOMASS, source], self.carbonEqu[YEAR_LATEST, EQ_BIOMASS, source])
                # clearBetaExcessPrev[source, target, :, :] = self.equation11and18(EQ_BIOMASS, source, TARGET)
                clearBetaExcessPrev[source, target] = equations.equation11and18(self.carbonExcess[YEAR_LATEST,EX_BIOMASS,source], self.transition_pft_fraction)
                # REFACTOR
                clearBetaExcessPrevSum[source, target, :] = self.equation12and19(clearBetaExcessPrev, source, target)
                clearBeta[source, target,:] = self.equation13and20(clearBetaEquPrev, clearBetaExcessPrevSum, source, target)
        
        # RUMUS 14 & RUMUS 15
        for i in range(len(targets)):
            clearBetaEquNext = self.equation14(sources, TARGETS[i])
            self.equation15(clearBetaEquNext, TARGETS[i], CLEARING)

        product_pools = [EX_PRODUCT1, EX_PRODUCT10, EX_PRODUCT100]
        product_clear_fracs = [CLEAR_P1, CLEAR_P10, CLEAR_P100]
        for i in range(len(product_pools)):
            for j in range(len(targets)):
                self.equation16_product(product_pools[i], TARGETS[j], clearBeta.sum(axis=0)[targets[j], :], product_clear_fracs[i])
        
        soil_pools = [EX_SOILRAPID, EX_SOILSLOW]
        for soil_pool in soil_pools:
            for j in range(len(TARGETS)):
                self.equation16_soil(soil_pool, TARGETS[j], clearBeta[targets[j],:])

    def soil_clearing(self, sources, targets, TARGETS):
        clearSigmaEquPrev, clearSigmaExcessPrev, clearSigmaExcessPrevSum, clearSigma = self.initializeSymbols(len(sources),len(targets), self.TRANS_NUM, self.PFT_NUM)
        SigmaEquNext = np.zeros((len(sources), len(targets), self.PFT_NUM))

        for source in sources:
            for j in range(len(targets)):
                clearSigmaEquPrev[source, targets[j], :] = self.equation10and17(EQ_BIOMASS, source, TARGETS[j])
                clearSigmaExcessPrev[source, targets[j], :, :] = self.equation11and18(EQ_BIOMASS, source, TARGETS[j])
                # REFACTOR
                clearSigmaExcessPrevSum[source, targets[j], :] = self.equation12and19(clearSigmaExcessPrev, source, targets[j])
                clearSigma[source, targets[j],:] = self.equation13and20(clearSigmaEquPrev, clearSigmaExcessPrevSum, source, targets[j])
        for source in sources:
            for j in range(len(TARGETS)):
                SigmaEquNext[source, targets[j],:] = self.equation21(EQ_SOILSLOW, source, TARGETS[j])
        soil_pools = [EX_SOILRAPID, EX_SOILSLOW]
        for soil_pool in soil_pools:
            for j in range(len(TARGETS)):
                self.equation22(soil_pool, TARGETS[j], targets[j], CLEARING, clearSigma, SigmaEquNext)

    def biomass_abandonment(self, sources, SOURCES, target, TARGET):
        # print()
        # print("abandon")
        """
        Menghitung perubahan jumlah karbon biomassa

        Parameters:
            sources (list(int)) : source cover type khusus untuk abandonment biomass
            SOURCES (list(int)) : source cover type absolut
            target (int) : target cover type khusus untuk abandonment biomass
            TARGET (int) : target cover type absolut
        """
        abandonBetaEquPrev, abandonBetaExcessPrev, abandonBetaExcessPrevSum, abandonBeta = self.initializeSymbols(len(sources), 1, self.TRANS_NUM, self.PFT_NUM)
        abandonBetaEquNext = np.zeros((self.PFT_NUM))

        for i in range(len(sources)):
            abandonBetaEquPrev[sources[i], target, :] = self.equation10and17(EQ_BIOMASS, SOURCES[i], TARGET)
            abandonBetaExcessPrev[sources[i], target, :, :] = self.equation11and18(EQ_BIOMASS, SOURCES[i], TARGET)
            abandonBetaExcessPrevSum[sources[i], target, :] = self.equation12and19(abandonBetaExcessPrev, sources[i], target)
            abandonBeta[sources[i], target,:] = self.equation13and20(abandonBetaEquPrev, abandonBetaExcessPrevSum, sources[i], target)
        abandonBetaEquNext = self.equation14(sources, TARGET)
        self.equation23(TARGET, ABANDON, abandonBetaEquNext, abandonBeta)

    def soil_abandonment(self, sources, SOURCES, target, TARGET):
        abandonSigmaEquPrev, abandonSigmaExcessPrev, abandonSigmaExcessPrevSum, abandonSigma = self.initializeSymbols(len(sources),1, self.TRANS_NUM, self.PFT_NUM)
        abandonSigmaEquNext = np.zeros((len(sources), self.PFT_NUM))

        for i in range(len(sources)):
            abandonSigmaEquPrev[sources[i], target, :] = self.equation10and17(EQ_BIOMASS, SOURCES[i], TARGET)
            abandonSigmaExcessPrev[sources[i], target, :, :] = self.equation11and18(EQ_BIOMASS, SOURCES[i], TARGET)
            # REFACTOR
            abandonSigmaExcessPrevSum[sources[i], target, :] = self.equation12and19(abandonSigmaExcessPrev, sources[i], target)
            abandonSigma[sources[i], target,:] = self.equation13and20(abandonSigmaEquPrev, abandonSigmaExcessPrevSum, sources[i], target)
        for i in range(len(sources)):
            abandonSigmaEquNext[sources[i],:] = self.equation21(EQ_SOILSLOW, SOURCES[i], TARGET)
        self.equation24(TARGET, ABANDON, abandonSigmaEquPrev, abandonSigmaEquNext.sum(axis=0))

    def biomass_others(self, source, SOURCE, target, TARGET):
        # print()
        # print("others")
        othersBetaEquPrev, othersBetaExcessPrev, othersBetaExcessPrevSum, othersBeta = self.initializeSymbols(1, 1, self.TRANS_NUM, self.PFT_NUM)
        othersBetaEquNext = np.zeros((self.PFT_NUM))

        othersBetaEquPrev[source, target, :] = self.equation10and17(EQ_BIOMASS, SOURCE, TARGET)
        othersBetaExcessPrev[source, target, :, :] = self.equation11and18(EQ_BIOMASS, SOURCE, TARGET)
        othersBetaExcessPrevSum[source, target, :] = self.equation12and19(othersBetaExcessPrev, source, target)
        othersBeta[source, target,:] = self.equation13and20(othersBetaEquPrev, othersBetaExcessPrevSum, source, target)
        othersBetaEquNext = self.equation14(SOURCE, TARGET)
        self.equation23(TARGET, OTHERS, othersBetaEquNext, othersBeta)

    def soil_others(self, source, SOURCE, target, TARGET):
        othersSigmaEquPrev, othersSigmaExcessPrev, othersSigmaExcessPrevSum, othersSigma = self.initializeSymbols(1,1, self.TRANS_NUM, self.PFT_NUM)
        othersSigmaEquNext = np.zeros((1, self.PFT_NUM))

        othersSigmaEquPrev[source, target, :] = self.equation10and17(EQ_BIOMASS, SOURCE, TARGET)
        othersSigmaExcessPrev[source, target, :, :] = self.equation11and18(EQ_BIOMASS, SOURCE, TARGET)
        # REFACTOR
        othersSigmaExcessPrevSum[source, target, :] = self.equation12and19(othersSigmaExcessPrev, source, target)
        othersSigma[source, target,:] = self.equation13and20(othersSigmaEquPrev, othersSigmaExcessPrevSum, source, target)
        othersSigmaEquNext[source,:] = self.equation21(EQ_SOILSLOW, SOURCE, TARGET)
        self.equation24(TARGET, ABANDON, othersSigmaEquPrev, othersSigmaEquNext.sum(axis=0))

    def biomass_harvest(self, SOURCES, target, TARGET):
        # print()
        # print('harvest')
        # harvestBetaEquPrev, harvestBetaExcessPrev, harvestBetaExcessPrevSum, harvestBeta = self.initializeSymbols(len(SOURCES), 1, self.TRANS_NUM, self.PFT_NUM)
        # harvestBetaEquNext = np.zeros((len(SOURCES), self.PFT_NUM))

        # for source in SOURCES:
        #     harvestBetaEquPrev[source, target, :] = self.equation10and17(EQ_BIOMASS, source, TARGET)
        #     harvestBetaExcessPrev[source, target, :, :] = self.equation25a(source, TARGET)
        #     harvestBetaExcessPrevSum[source, target, :] = self.equation12and19(harvestBetaExcessPrev, source, target)
        # harvestBetaEquNext = self.equation14(SOURCES, TARGET)
        # self.equation26(harvestBetaEquNext, harvestBetaExcessPrevSum)
        # harvestBeta = self.equation27(harvestBetaEquPrev, harvestBetaExcessPrevSum)
        # self.equation28(harvestBeta.sum(axis=1))
        # self.equation29(harvestBeta.sum(axis=1))
        BetaHarvestExcessPrev = np.zeros((len(SOURCES), self.TRANS_NUM, self.PFT_NUM))

        # EQUATION 10
        BetaHarvestEquPrev = self.area * self.harvest_transition_pft_fraction * c.carbonEquDense[EQ_BIOMASS,VIRGIN:SECOND+1]
        self.carbonEqu[YEAR_LATEST,EQ_BIOMASS,VIRGIN:SECOND+1,:] = self.carbonEqu[YEAR_LATEST,EQ_BIOMASS,VIRGIN:SECOND+1,:] - BetaHarvestEquPrev
        
        for source in SOURCES:
            # RUMUS 25a
            BetaHarvestExcessPrev[source] = self.carbonExcess[YEAR_LATEST, EX_BIOMASS,source, :,:] * np.tile(self.harvest_transition_pft_fraction[source, np.newaxis,:], (1,self.TRANS_NUM, 1)) / self.temp_harvest_cover_fraction[source]
        # RUMUS 25b
        BetaHarvestExcessPrevSum = BetaHarvestExcessPrev.sum(axis=1)

        # RUMUS 14 harvest
        BetaHarvestEquNext = self.area * self.harvest_transition_pft_fraction.sum(axis=0) * c.carbonEquDense[EQ_BIOMASS, SECOND]
        self.carbonEqu[YEAR_LATEST,EQ_BIOMASS, SECOND, :] = self.carbonEqu[YEAR_LATEST,EQ_BIOMASS, SECOND, :] + BetaHarvestEquNext

        # RUMUS 26
        self.carbonExcess[YEAR_LATEST, EX_BIOMASS, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, EX_BIOMASS,SECOND,HARVEST,:] - (BetaHarvestEquNext + BetaHarvestExcessPrevSum).sum(axis=0)

        # RUMUS 27
        BetaHarvest = BetaHarvestEquPrev + BetaHarvestExcessPrevSum

        # RUMUS 28
        self.carbonExcess[YEAR_LATEST, EX_SOILRAPID, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, EX_SOILRAPID, SECOND, HARVEST, :] + (BetaHarvest * c.harvest[HARVEST_SOIL_VIRGIN: HARVEST_SOIL_SECOND+1]).sum(axis=0)

        # RUMUS 29
        remainder = 1 - c.harvest[HARVEST_SOIL_VIRGIN: HARVEST_SOIL_SECOND+1]
        products = [EX_PRODUCT1, EX_PRODUCT10, EX_PRODUCT100]
        product_constants = [HARVEST_P1, HARVEST_P10, HARVEST_P100]
        for i in range(len(products)):
            self.carbonExcess[YEAR_LATEST, products[i], SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, products[i], SECOND, HARVEST, :] + (BetaHarvest * remainder).sum(axis=0) * c.harvest[product_constants[i]]
        
    def soil_harvest(self, SOURCES, target, TARGET):
        # harvestSigmaEquPrev, harvestSigmaExcessPrev, harvestSigmaExcessPrevSum, harvestSigma = self.initializeSymbols(len(SOURCES), 1, self.TRANS_NUM, self.PFT_NUM)
        # harvestSigmaEquNext = np.zeros((len(SOURCES), self.PFT_NUM))

        # for source in SOURCES:
        #     harvestSigmaEquPrev[source, target, :] = self.equation10and17(EQ_SOILSLOW, source, TARGET)
        #     harvestSigmaExcessPrev[source, target,:] = self.equation30a(source)
        # harvestSigmaExcessPrevSum = harvestSigmaExcessPrev.sum(axis=2)
        # harvestSigma = harvestSigmaEquPrev + harvestSigmaExcessPrevSum

        # RUMUS 17 Harvest
        SigmaHarvestEquPrev = self.area * self.harvest_transition_pft_fraction * c.carbonEquDense[EQ_SOILSLOW,VIRGIN:SECOND+1]
        self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW,VIRGIN:SECOND+1,:] = self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW,VIRGIN:SECOND+1,:] - SigmaHarvestEquPrev

        # RUMUS 30
        SigmaHarvestExcessPrev = np.zeros((len(SOURCES), self.TRANS_NUM, self.PFT_NUM))
        for source in SOURCES:
            SigmaHarvestExcessPrev[source] = self.carbonExcess[YEAR_LATEST, EX_BIOMASS,source, :,:] * np.tile(self.harvest_transition_pft_fraction[source, np.newaxis,:], (1,self.TRANS_NUM, 1)) / self.temp_harvest_cover_fraction[source]
        SigmaHarvestExcessPrevSum = SigmaHarvestExcessPrev.sum(axis=1)

        # RUMUS 31
        SigmaHarvest = SigmaHarvestEquPrev + SigmaHarvestExcessPrevSum
        # RUMUS 21 Harvest
        SigmaHarvestEquNext = self.area * self.harvest_transition_pft_fraction.sum(axis=0) * c.carbonEquDense[EQ_SOILSLOW, SECOND]
        self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW, SECOND, :] = self.carbonEqu[YEAR_LATEST,EQ_SOILSLOW, SECOND, :] + SigmaHarvestEquNext
        # RUMUS 32
        SigmaRapid = np.maximum(0, SigmaHarvest.sum(axis=0) - self.area * self.harvest_transition_pft_fraction.sum(axis=0) * c.harvestMinimumSoil)
        self.carbonExcess[YEAR_LATEST, EX_SOILRAPID, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, EX_SOILRAPID, SECOND, HARVEST, :] + SigmaRapid
        # RUMUS 33
        self.carbonExcess[YEAR_LATEST, EX_SOILSLOW, SECOND, HARVEST, :] = self.carbonExcess[YEAR_LATEST, EX_SOILRAPID, SECOND, HARVEST, :] - SigmaRapid + SigmaHarvestEquPrev.sum(axis=0) + SigmaHarvestEquNext

    def initializeSymbols(self, len_sources, len_targets, len_TRANS, len_PFT):
        symbol_without_TRANS = np.zeros((len_sources, len_targets, len_PFT))
        symbol_with_TRANS = np.zeros((len_sources, len_targets, len_TRANS, len_PFT))
        return (symbol_without_TRANS,
                symbol_with_TRANS,
                symbol_without_TRANS.copy(),
                symbol_without_TRANS.copy())

    def createHarvestCoverDivisor(self):
        divisor = self.cover_fraction.copy()
        divisor[divisor==0] = 1 
        divisor = np.tile(divisor[VIRGIN:SECOND+1,np.newaxis],(1,11))
        return divisor

    def createTransitionAxis(self):
        
        
    #     result = np.zeros((2,2,self.FOREST_IDX))
        
    #     result[VIRGIN, FORESTs,:] = self.cover_pft_fraction[VIRGIN, 0:self.FOREST_IDX]
    #     result[VIRGIN, NONFOREST,:] = self.cover_pft_fraction[VIRGIN, self.FOREST_IDX:]
    #     # result[VIRGIN, NONFOREST, :-1] = self.cover_pft_fraction[VIRGIN, self.FOREST_IDX:]
    #     result[SECOND, FOREST,:] = self.cover_pft_fraction[SECOND, 0:self.FOREST_IDX]
    #     result[SECOND, NONFOREST,:] = self.cover_pft_fraction[SECOND, self.FOREST_IDX:]
    #     # result[SECOND, NONFOREST,:-1] = self.cover_pft_fraction[SECOND, self.FOREST_IDX:]
    #     return result
        
    
        
    
        
    
        
