# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 11:04:13 2023

@author: Pierre
"""


import numpy as np
import xarray as xr
import warnings
import extract_mat
import time
import os
# import rioxarray


# No warnings
warning_off = True
if warning_off:
    warnings.filterwarnings("ignore", category=RuntimeWarning)
np.set_printoptions(suppress=True)  # Put off scientific notation

now = time.time()

# Correct paths and options
platform = "local"
if platform=="local":
    path_data = '/home/zeigerp/Documents/Postdoc_IGE/data'
    name = 'Antarctica'
    grid = xr.open_dataset(path_data+'/SMOS_data/%s_grid_12km_lonlat.nc'%name)
    # Reconstruct
    print("Starting reconstruction")
    ds_TB = extract_mat.read_mat_day([path_data + '/SMOS_data/mat', grid, 'West'])
    
elif platform=="dahu":
    path_data = '/bettik/PROJECTS/pr-snowem/zeigerp'
    name = os.environ["NAME"]
    subset = os.environ["SUBSET"]
    assert (name=='Greenland') | (name=='Antarctica'), "$NAME variable should be one of: 'Greenland', 'Antarctica'"
    grid = xr.open_dataset(path_data+'/mask/%s_grid_12km_lonlat.nc'%name)
    # Reconstruct
    print("Starting reconstruction")
    ds_TB = extract_mat.reconstruct_from_mat(path_data, grid, subset=None) 

print('Finished. Total time ellapsed: %.2f hours'%((time.time()-now)/3600))





