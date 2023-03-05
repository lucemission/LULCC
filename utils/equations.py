import numpy as np
from index import YEAR_LATEST

def equation7(carbon_equilibrium, area, carbon_equilibrium_density):
    return carbon_equilibrium / (area * carbon_equilibrium_density)

def equation9(cover_pft_fraction, transition_fraction, cover_fraction):
    return cover_pft_fraction * transition_fraction / cover_fraction

def equation10and17(area, transition_pft_fraction, carbon_equilibrium_density, carbon_equilibrium):
    emitted_carbon_equilibrium = area * transition_pft_fraction[source, target, :] * carbon_equilibrium_density[pool, source,:]
    carbon_equilibrium[YEAR_LATEST,pool, source, :] = carbon_equilibrium[YEAR_LATEST,pool, source, :] - emitted_carbon_equilibrium
    return emitted_carbon_equilibrium

def equation11and18(carbon_excess, transition_pft_fraction, cover_pft_fraction):
    emitted_carbon_excess = carbon_excess * transition_pft_fraction / cover_pft_fraction
    carbon_excess = carbon_excess - emitted_carbon_excess
    return emitted_carbon_excess
