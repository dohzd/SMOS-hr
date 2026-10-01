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

from edit_nc import edit_tb_attributes

name = 'Antartica'
export_mode = 'yearly'
incidence = 40
path = '/bettik/PROJECTS/pr-snowem/zeigerp/SMOS_data/rSIR-enhanced/single_incidence_%s/%s'%(incidence, name)   

if name == "Antartica":
    # melt4d = xr.open_zarr("https://snow.univ-grenoble-alpes.fr/opendata/melt-4D-Antarctic-12km.zarr")
    melt4d = xr.open_dataset('/bettik/PROJECTS/pr-snowem/zeigerp/SMOS_data/melt4d_lon_lat.nc')
else:
    melt4d = xr.open_dataset('/bettik/PROJECTS/pr-snowem/zeigerp/SMOS_data/greenland_lon_lat.nc')

if export_mode == 'yearly':
    
    years = glob.glob(path + '/20*')
    years = np.sort([int(y.split('/')[-1]) for y in years])
    
    for i, y in enumerate(years):
        print('%s...'%y)
        
        files = np.sort(glob.glob(path + '/%s/TB_rSIR_%s_*.nc'%(y, name)))
        times = np.asarray([datetime.strptime(f.split('_')[-4].split('.')[0],'%Y-%m-%d') for f in files])
        xr_all = xr.open_mfdataset(files, combine = 'nested', concat_dim = 'time')#.chunk(time = 366)
        xr_all['time'] = times
        xr_sel = xr_all[['TB_H_morning', 'TB_H_afternoon', 'TB_V_morning', 'TB_V_afternoon']].sel(iterations = 10).astype('float32').compute()

        # Change attributes dtype, names and format the dataset
        xr_sel = edit_tb_attributes(xr_sel, name)     
        
        # Add lon and lat of the center of each pixel from melt4d files
        xr_sel['lon'] = melt4d['lon']
        xr_sel['lat'] = melt4d['lat']
    
        # Compression
        comp = dict(zlib=True, complevel=5)    
        encoding = {var: comp for var in xr_sel.data_vars}
        
        if name == "Antartica":
            xr_sel.to_netcdf(path + '/SMOS-rSIR_enhanced-%s-Antarctica-%s.nc'%(incidence, y), encoding=encoding)     
        else:   
            xr_sel.to_netcdf(path + '/SMOS-rSIR_enhanced-%s-%s-%s.nc'%(incidence, name, y), encoding=encoding)        
        # xr_sel.chunk(time = -1, x = 100, y = 100).to_zarr(path + '/%s/SMOS_rSIR-enhanced-%s_%s_%s.zarr'%(y, incidence, name, y))
        
elif export_mode == 'all':
    
    files = np.sort(glob.glob(path + '/20*/TB_rSIR_%s_*.nc'%name))
    times = np.asarray([datetime.strptime(f.split('_')[-4].split('.')[0],'%Y-%m-%d') for f in files])
    
    xr_all = xr.open_mfdataset(files, combine = 'nested', concat_dim = 'time').chunk(dict(time = -1, x = 30, y = 30))
    xr_all['time'] = times
    xr_sel = xr_all[['TB_H_morning', 'TB_H_afternoon', 'TB_V_morning', 'TB_V_afternoon']].sel(iterations = 10).astype('float32')#.compute()
    
    # Change attributes dtype, names and format the dataset
    xr_sel = edit_tb_attributes(xr_sel, name)

    # Add lon and lat of the center of each pixel from melt4d files
    xr_sel['lon'] = melt4d['lon']
    xr_sel['lat'] = melt4d['lat']
    
    # Compression
    comp = dict(zlib=True, complevel=5)    
    encoding = {var: comp for var in xr_sel.data_vars}
            
    if name == 'Antartica':
        name = 'Antarctica'
    xr_sel.to_netcdf(path + '/SMOS-rSIR_enhanced-%s-%s-2010_2024.nc'%(incidence, name), encoding=encoding)        










