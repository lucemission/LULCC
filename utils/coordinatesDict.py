# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 23:08:28 2023

@author: Khalid
"""
import os
import re
import netCDF4 as nc
import numpy as np
from index import VIRGIN, SECOND, FORESTs, F_SBH, F_SBH2, F_SBH3, F_VBH, F_VBH2, NATURAL_PFTs

from Coordinate import Coordinate

def coordinatesDict(data_path, start_year, duration, input_coordinate):
    """
    Membuat dictionary berisi kelas Coordinate untuk tiap pasangan koordinat latitude dan longitude
    
    Parameters:
        data_path (string): alamat data
        start_year (int) : tahun dimulainya bookkeeping
        duration (int) : durasi berjalannya bookkeeping
        input_coordinate (tuple(int)): koordinat tertentu yang ingin dilakukan bookkeeping
        
    Returns:
        coordinates (dict(Coordinate)) : dictionary berisi kelas Coordinate
    """
    coordinates = dict()
    
    lat_coordinates = np.linspace(89.75, -89.75, 360)
    lon_coordinates = np.linspace(-179.75, 179.75, 720)
    
    states = statesReader(data_path, start_year, duration)
    pft = nc.Dataset(f"{data_path}/PFT/LAND_COV_scen_all.nc")['vegfract'][:duration+1, :NATURAL_PFTs].data
    # pft_1500 = pft_1500_nc['vegfract'][0,:11,180,204].data
    # print(pft_1500.shape)
    
    input_lat, input_lon = input_coordinate
    filename = f"lat{input_lat}lon{input_lon}.lu"
    
    # for filename in os.listdir(f"{data_path}/lu"):  
    # transitions = np.loadtxt(f"{data_path}/lu/{filename}", skiprows=1, usecols=tuple(np.arange(1,12)))
    transition_columns = np.append(np.arange(1,12), np.array([13,15,17,19,21]))
    transition = np.loadtxt(f"{data_path}/lu/{filename}", skiprows=1, usecols=transition_columns)
    # print(filename)
    [lat, lon] = extractLatLon(filename)
    
    lat_idx = np.where(lat_coordinates == lat)
    lon_idx = np.where(lon_coordinates == lon)
                     
    string = re.sub("(\.|lu)", "",filename)
    coordinates[f"{re.sub('-', '_', string)}"] = Coordinate(
                                                    transition=cubeTransitions(transition),
                                                    harvest_transition = cubeHarvestTransition(transition),
                                                    initial_states = extract(states, lat_idx, lon_idx),
                                                    pft = extract(pft, lat_idx, lon_idx),
                                                    duration=duration)
    return coordinates

def extractLatLon(filename):
    lat = re.search('lat(.+?)lon', filename).group(1)
    lon = re.search('lon(.+?)\.lu', filename).group(1)
    return (float(lat), float(lon))

def extract(fractions, lat_idx, lon_idx):
    # return np.squeeze(initial_states[:,lat_idx, lon_idx])
    return sumToOne(np.squeeze(fractions[:,:,lat_idx, lon_idx]))

def sumToOne(fractions):
    return fractions/fractions.sum(axis=1)[:, np.newaxis]

def statesReader(data_path, start_year, duration):
    """
    Mengambil data cover (a.k.a states) fractions
    
    Parameters:
        data_path (string): alamat data
        start_year (int): tahun dimulainya bookkeeping
        duration (int): durasi berjalannya bookkeeping
        
    Returns:
        states (array(int)): pecahan luas cover tiap tahunnya
    """
    
    states_path    = f"{data_path}/updated_states"
    simulated_years = start_year + np.linspace(0, duration, duration+1, dtype=int)
    for year in simulated_years:
        if year == start_year:
            initial_virgin = np.loadtxt(f'{states_path}/gothr.{year}.txt')
            initial_sec    = np.loadtxt(f'{states_path}/gsecd.{year}.txt')
            initial_pasture= np.loadtxt(f'{states_path}/gpast.{year}.txt')
            initial_crop   = np.loadtxt(f'{states_path}/gcrop.{year}.txt')
            states = np.array([[initial_virgin, initial_sec, initial_pasture, initial_crop]])
        else:
            virgin = np.loadtxt(f'{states_path}/gothr.{year}.txt')
            sec    = np.loadtxt(f'{states_path}/gsecd.{year}.txt')
            pasture= np.loadtxt(f'{states_path}/gpast.{year}.txt')
            crop   = np.loadtxt(f'{states_path}/gcrop.{year}.txt')
            append_states = np.array([[virgin, sec, pasture, crop]])
            states = np.append(states,append_states, axis=0)
    return states

def cubeTransitions(transitions1D):
    transitions2D = np.zeros((transitions1D.shape[0],4,4))
    # transitions2D[:,0,1] = transitions1D[:,10]
    transitions2D[:,0,2] = transitions1D[:,3]
    transitions2D[:,0,3] = transitions1D[:,4]
    transitions2D[:,1,2] = transitions1D[:,8]
    transitions2D[:,1,3] = transitions1D[:,6]
    transitions2D[:,2,0] = transitions1D[:,2]
    transitions2D[:,2,1] = transitions1D[:,9]
    transitions2D[:,2,3] = transitions1D[:,1]
    transitions2D[:,3,0] = transitions1D[:,5]
    transitions2D[:,3,1] = transitions1D[:,7]
    transitions2D[:,3,2] = transitions1D[:,0]
    return transitions2D

def cubeHarvestTransition(transitions1D):
    harvestTransition2D = np.zeros((transitions1D.shape[0],2,11))
    harvestTransition2D[:, VIRGIN, :FORESTs] = np.tile(transitions1D[:,F_VBH,np.newaxis], (1,6))
    harvestTransition2D[:, VIRGIN, FORESTs: ] =  np.tile(transitions1D[:,F_VBH2, np.newaxis], (1,5))
    harvestTransition2D[:, SECOND, :FORESTs] =  np.tile(sumSBH(transitions1D)[:,np.newaxis], (1,6))
    harvestTransition2D[:, VIRGIN, FORESTs:] =  np.tile(transitions1D[:,F_SBH3,np.newaxis], (1,5))
    return harvestTransition2D

def sumSBH(transitions1D):
    return transitions1D[:,F_SBH]+transitions1D[:,F_SBH2]

def simulation():
    transition = np.array([[[0,0.2,0.1,0.05],
                       [0,0,0.08,0.1],
                       [0,0.1,0,0.03],
                       [0,0.07,0.01,0]]])
    harvest_transition = np.array([[[0.15,0.05],
                                   [0.12,0.08]]])
    cover = np.array([0.3,0.5,0.1,0.1])
    pft = np.array([0.8,0.2])
    return Coordinate(transition, harvest_transition, cover, pft)

# def pftReader():
#     return
    
    