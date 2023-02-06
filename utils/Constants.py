# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 01:19:46 2023

@author: Khalid
"""
import numpy as np

carbonEquDense = np.array([
      #PFT 1     2     3     4     5     6  7  8   9  10  11
    [# Biomass
      [  200,  160,  160,  135,   90,   90,27,27,  7, 18,  3], # Primary   - Veg
      [  150,  120,  120,  100,   68,   68,27,27,  7, 18,  0], # Secondary - Veg
      [   18,   18,    7,    7,    7,    7,18, 7,  7, 18,  7], # Pasture   - Veg
      [    5,    5,    5,    5,    5,    5, 5, 5,  3,  5,  1]], # Crop      - Veg
    [# Soil
      [  117,  117,  134,  134,  206,  206,69,69,189, 42,204], # Primary   - Soil
      [   88,   88,  120,  120,  185,  185,69,69,189, 42,204], # Secondary - Soil
      [87.75,87.75,100.5,100.5,154.5,154.5,69,69,189, 42,204], # Pasture   - Soil
      [   58,   58,   67,   67,  103,  103,34,34, 94, 21,101]]  # Crop      - Soil
])

# carbonEquDense = np.array([
#       #PFT 1     2     3     4     5     6  7  8   9  10  11
#     [# Biomass
#       [  200,  160], # Primary   - Veg
#       [  150,  120], # Secondary - Veg
#       [   18,   18], # Pasture   - Veg
#       [    5,    5]], # Crop      - Veg
#     [# Soil
#       [  117,  117], # Primary   - Soil
#       [   88,   88], # Secondary - Soil
#       [87.75,87.75], # Pasture   - Soil
#       [   58,   58]]  # Crop      - Soil
# ])

# harvestEquDense = np.array([
#     [# Biomass
#       [  [200],  [160]], # Primary   - Veg
#       [  [150],  [120]] # Secondary - Veg
#       ],
#     [# Soil
#       [  [117],  [117]], # Primary   - Soil
#       [   [88],   [88]]] # Secondary - Soil
#     ])

harvest = np.array([
    #PFT 1    2    3    4    5    6    7    8    9   10   11
    [ 0.9, 0.9, 0.4, 0.4, 0.4, 0.4,   1,   1,   1,   1,   1], # Fraction of roundwood assigned to decay pools after harvest, 1yr
    [0.04,0.04,0.24,0.24,0.24,0.24,   0,   0,   0,   0,   0], # Fraction of roundwoodassigned to decay pools after harvest, 10yr
    [0.06,0.06,0.36,0.36,0.36,0.36,   0,   0,   0,   0,   0], # Fraction of roundwood assigned to decay pools after harvest, 100yr
    [0.79,0.86,0.81,0.78,0.87,0.87,0.86,0.78,0.78,0.86,0.87], # Fraction of vegetation carbon transferred dead to soil at clearing for primary forest
    [0.71,0.81,0.75,0.70,0.82,0.82,0.81,0.70,0.70,0.81,0.82], # Fraction of vegetation carbon transferred dead to soil at clearing for secondary forest
    [  76,  76,  67,  67, 175, 175,  45,  35,  95,  27, 173] # Minimum soil C following harvest
#     # [   5,   5,  10,  10,  15,  15,   5,  10,  10,   5,  15] # 
])

timeMinimumHarvest = np.array([   5,   5,  10,  10,  15,  15,   5,  10,  10,   5,  15])  # Time of soil carbon to reach minimum (as found following harvest)

clearing = np.array([
    #PFT 1    2    3    4    5    6   7   8   9  10  11
    [ 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,0.4,0.4,0.5,0.5,0.4], # Fraction of vegetation carbon assigned to decay pools after clearing, 1yr
    [0.27,0.27, 0.2, 0.2, 0.2, 0.2,0.1,0.1,  0,  0,0.1], # Fraction of vegetation carbon assigned to decay pools after clearing, 10yr
    [   0,   0,0.07,0.07,0.07,0.07,  0,  0,  0,  0,  0], # Fraction of vegetation carbon assigned to decay pools after clearing, 100yr
    [0.33,0.33,0.33,0.33,0.33,0.33,0.5,0.5,0.5,0.5,0.5], # Fraction of vegetation carbon transferred dead to soil at clearing (1-SUM(product pools))
    [ 0.8, 0.8,0.81, 0.81,0.45,0.45,0.8,0.8,0.8,0.81,0.81],
    [ 0.6, 0.6,0.75,0.75, 0.3, 0.3,0.8,0.8,0.8,0.81,0.81]
    # Soil carbon after initial, rapid loss after clearing
])

# clearing = np.array([
#     #PFT 1    2    3    4    5    6   7   8   9  10  11
#     [ 0.4, 0.3], # Fraction of vegetation carbon assigned to decay pools after clearing, 1yr
#     [ 0.2, 0.3], # Fraction of vegetation carbon assigned to decay pools after clearing, 10yr
#     [0.07,0.05], # Fraction of vegetation carbon assigned to decay pools after clearing, 100yr
#     [0.33,0.35], # Fraction of vegetation carbon transferred dead to soil at clearing (1-SUM(product pools))
#     [0.81, 0.81],
#     [0.75,0.75]
#     # Soil carbon after initial, rapid loss after clearing
# ])

# recovery = np.array([
#  #PFT 1  2  3  4  5  6  7  8  9 10 11
#     [50,50,50,50,50,50,25,50,10, 5,50], # Time required for biomass carbon to recover (to secondary land) after abandoned
#     [15,15,40,40,35,35,15,45,45,15,45] # Time required for soil carbon to recover (to secondary land) after abandoned
# ])


# Timescale = np.array([
# #    #7,4,4,11
# # BIOMASS
#     [ # virgin
#      [# clearing
#       [0,0,0,0,0,0,0,0,0,0,0],
#       # harvest
#       [0,0,0,0,0,0,0,0,0,0,0],
#       # abandonment
#       [0,0,0,0,0,0,0,0,0,0,0],
#       # others
#       [0,0,0,0,0,0,0,0,0,0,0]],
#      # secondary
#      [# clearing
#       [0,0,0,0,0,0,0,0,0,0,0],
#       # harvest
#       [26.7,26.7,26.7,26.7,26.7,26.7,13.35,26.7,5.34,2.67,26.7],
#       # abandonment
#       [26.7,26.7,26.7,26.7,26.7,26.7,13.35,26.7,5.34,2.67,26.7],
#       # others
#       [0,0,0,0,0,0,0,0,0,0,0]],
#      # pasture
#      [# clearing
#       [0,0,0,0,0,0,0,0,0,0,0],
#       # harvest
#       [0,0,0,0,0,0,0,0,0,0,0],
#       # abandonment
#       [0,0,0,0,0,0,0,0,0,0,0],
#       # others
#       [0,0,0,0,0,0,0,0,0,0,0]],
#      # crop
#      [# clearing
#       [0,0,0,0,0,0,0,0,0,0,0],
#       # harvest
#       [0,0,0,0,0,0,0,0,0,0,0],
#       # abandonment
#       [0,0,0,0,0,0,0,0,0,0,0],
#       # others
#       [0,0,0,0,0,0,0,0,0,0,0]]],
# # SS
#     [# virgin
#      [ # clearing
#       [],
#       # harvest
#       [],
#       # abandonment
#       [],
#       # others
#       []],
#      # secondary
#      [# clearing
#       [],
#       # harvest
#       [],
#       # abandonment
#       [],
#       # others
#       []],
#      # pasture
#      [# clearing
#       [],
#       # harvest
#       [],
#       # abandonment
#       [],
#       # others
#       [1.602, 1.602, 8.01, 8.01, 26.7, 26.7, 1.602, 8.01, 8.01, 1.602, 8.01]],
#      # crop
#      [# clearing
#       [],
#       # harvest
#       [],
#       # abandonment
#       [],
#       # others
#       [1.602, 1.602, 8.01, 8.01, 26.7, 26.7, 1.602, 8.01, 8.01, 1.602, 8.01]]
#      ],
# # SR
#     [# virgin
#      [ # clearing
#       [],
#       # harvest
#       [],
#       # abandonment
#       [],
#       # others
#       []],
#      # secondary
#      [# clearing
#       [],
#       # harvest
#       [],
#       # abandonment
#       [],
#       # others
#       []],
#      # pasture
#      [# clearing
#       [1.602, 1.602, 8.01, 8.01, 26.7, 26.7, 1.602, 8.01, 8.01, 1.602, 8.01],
#       # harvest
#       [],
#       # abandonment
#       [],
#       # others
#       []],
#      # crop
#      [# clearing
#       [1.602, 1.602, 8.01, 8.01, 26.7, 26.7, 1.602, 8.01, 8.01, 1.602, 8.01],
#       # harvest
#       [],
#       # abandonment
#       [],
#       # others
#       []]
#      ],
# # P1 
#     [],
# # P10
#     [],
# # P100
#     []

    
# #   virgin
# #       clearing
    
# #       harvest
# #       abandonment
# #       others     
  
# ### secondary
# ### pasture
# ### crop
#     [],
# ])