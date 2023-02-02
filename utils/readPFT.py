# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import netCDF4 as nc
data = nc.Dataset('../data/PFT/1700-dst/LAND_COV_scen_1700_1992_293.nc')
pft1 = data['vegfract'][:,0,:,:].data
pft0mask = data['vegfract'][:,0,:,:].mask
#print(PFT['vegtype'][:,0,:,:].data)
pft2 = data['vegfract'][:,1,:,:].data
pft3 = data['vegfract'][:,2,:,:].data

