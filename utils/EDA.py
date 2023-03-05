import numpy as np
import pandas as pd
import netCDF4 as nc

geo_path = "../../data/staticData_quarterdeg.nc"
geo_nc = nc.Dataset(geo_path)
