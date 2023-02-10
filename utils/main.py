# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 11:24:03 2023

@author: Khalid
"""
# from Coordinate import Coordinate
from coordinatesDict import coordinatesDict, simulation
# import os

data_04_path = "../../data/1500-2004"
input_coordinate = (50.25, 10.25)
start_year = 1500
duration = 2
coordinates = coordinatesDict(data_04_path, start_year, duration, input_coordinate)

coordinates['lat5025lon1025'].run()
# coordinates['lat5025lon1025'].clearing()

# print(coordinates['lat_025lon_4825'].transitions[0,0])
# simulation = simulation()
# simulation.updateFractions()
# simulation.clearing()