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

path_data = '/bettik/PROJECTS/pr-snowem/zeigerp'
#name = 'Greenland'
name = os.environ["NAME"]

# Read AMSR and define 12,5 km grid
if name == 'Antartica':
    amsr = xr.open_dataset(path_data+'/AMSR/melt-2002-2023.nc')
    grid = amsr.isel(time=6502).snow_status_wet_dry_19H_ASC_filter
elif name == 'Greenland':
    grid = xr.open_dataset(path_data + '/mask/grid_AMSR_Greenland.nc')
    #grid = rioxarray.open_rasterio(path_data + '/mask/mask_glaciated_area_greenland_12km_bool.tif')


print("Starting reconstruction")

#%% Read L1C data and reconstuct TB images
now = time.time()
ds_TB = extract_mat.reconstruct_from_mat(path_data, grid, subset=None) 
print('Total time ellapsed: %.2f hours'%((time.time()-now)/3600))





