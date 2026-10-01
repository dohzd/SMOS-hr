#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:58:17 2024

@author: zeigerp
"""

import numpy as np
import xarray as xr
import glob
from datetime import datetime

# from detect_melt import compute_contamination_mask
from edit_nc import edit_melt_attributes


# Initially from detect_melt
def compute_contamination_mask(melt):    # Years should be a sequence containing all the m
    years = np.unique(melt.time.values.astype('datetime64[Y]').astype(int) + 1970)
    winter_melt = 0
    ndays = 0
    for y in years[:-1]:
        m = (melt > 0).sel(time = slice(f'{y}-06-15', f'{y}-08-15')).sum(['time'])
        winter_melt = winter_melt + m   
        ndays = ndays + len(melt.sel(time = slice(f'{y}-06-15', f'{y}-08-15')).time)  
    return winter_melt / ndays, years


name = 'Greenland'
path = '/bettik/PROJECTS/pr-snowem/zeigerp/SMOS_data/melt/%s/'%name
    
files = np.sort(glob.glob(path + 'SMOS_melt_12km_%s*.nc'%name))
    
ds_melt = xr.open_mfdataset(files, combine = 'nested', concat_dim = 'time')#.chunk(dict(time = -1, x = 30, y = 30))

# Compute contamination mask
ds_melt = ds_melt.sel(time = slice('2010-04-01', '2024-03-31'))
mask_contamination, years = compute_contamination_mask(ds_melt['snow_status_wet_dry_smos_daily'])
mask_contamination = mask_contamination.compute()

# Output variables and apply masks
if name == 'Antartica': 	# Contamination mask does not apply to Greenland
    ds_melt['snow_status_wet_dry_smos_morning'] = ds_melt['snow_status_wet_dry_smos_morning'].where(mask_contamination <= 5 / 100, other = -10)
    ds_melt['snow_status_wet_dry_smos_afternoon'] = ds_melt['snow_status_wet_dry_smos_afternoon'].where(mask_contamination <= 5 / 100, other = -10)
ds_melt['snow_status_wet_dry_smos_morning'] = ds_melt['snow_status_wet_dry_smos_morning'].where(ds_melt.mask_TBV_daily, other = -10)
ds_melt['snow_status_wet_dry_smos_afternoon'] = ds_melt['snow_status_wet_dry_smos_afternoon'].where(ds_melt.mask_TBV_daily, other = -10)

xr_sel = ds_melt[['snow_status_wet_dry_smos_morning', 'snow_status_wet_dry_smos_afternoon']]

# Add lon and lat of the center of each pixel from melt4d files
if name == "Antartica":
    melt4d = xr.open_dataset('/bettik/PROJECTS/pr-snowem/zeigerp/SMOS_data/melt4d_lon_lat.nc')
else:
    melt4d = xr.open_dataset('/bettik/PROJECTS/pr-snowem/zeigerp/SMOS_data/greenland_lon_lat.nc')
xr_sel['lon'] = melt4d['lon']
xr_sel['lat'] = melt4d['lat']
          
# Add attributes to the dataset
xr_sel = edit_melt_attributes(xr_sel, name)

# Compression
comp = dict(zlib=True, complevel=5)    
encoding = {var: comp for var in xr_sel.data_vars}
      
if name == 'Antartica':
    name = 'Antarctica'
xr_sel.to_netcdf(path + '/SMOS-melt-12km-%s-2010_2024.nc'%name, encoding=encoding)        



