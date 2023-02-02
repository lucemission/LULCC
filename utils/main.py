# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 11:24:03 2023

@author: Khalid
"""
# from Coordinate import Coordinate
from coordinatesDict import coordinatesDict, simulation
# import os

data_04_path = "../../latihan"
coordinates = coordinatesDict(data_04_path)

coordinates['lat_025lon_7775'].updateFractions()
coordinates['lat_025lon_7775'].clearing()

# print(coordinates['lat_025lon_4825'].transitions[0,0])
# simulation = simulation()
# simulation.updateFractions()
# simulation.clearing()