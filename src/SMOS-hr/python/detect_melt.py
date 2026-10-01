# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 11:28:30 2023

@author: Pierre
"""


import numpy as np
import xarray as xr
import glob
from datetime import datetime
import time
import os
import rioxarray

import melt4d.base_algorithm as algo
from melt4d.constants import WMO_MELT_VAR


def smos_melt_loop(path_data, name = 'Antartica'):

    # Open AMSR melt mask (<1700 m elevation)
    if name == 'Antartica':
        amsr_mask = ~np.isnan(xr.open_dataset(path_data+'/AMSR/melt-2002-2023.nc').snow_status_wet_dry_19H_ASC_filter.mean('time')) 	#.chunk(dict(time=1000, x=100, y=100))
        
    elif name == 'Greenland':
        amsr_mask = rioxarray.open_rasterio(path_data+'/mask/GrIS_fraction_12km_80%_clean.tif')
        
    # All melt years
    years = glob.glob(path_data + '/SMOS_data/rSIR-enhanced/single_incidence_40/%s/20*/'%name)
    years = np.sort([int(y.split('/')[-2]) for y in years])
    
    for i, y in enumerate(years[:-1]):      # Loop over the melt years
    
        # 1-year time series from April 01 to March 31
        files1 = np.sort(glob.glob(path_data + '/SMOS_data/rSIR-enhanced/single_incidence_40/%s/%s/TB_rSIR_%s_*.nc'%(name, str(y), name)))
        files2 = np.sort(glob.glob(path_data + '/SMOS_data/rSIR-enhanced/single_incidence_40/%s/%s/TB_rSIR_%s_*.nc'%(name, str(y + 1), name)))
        files = np.concatenate((files1, files2))
        times = np.asarray([datetime.strptime(f.split('_')[-4].split('.')[0],'%Y-%m-%d') for f in files])
        time_mask = (times < datetime(y + 1, 4, 1)) & ((times > datetime(y, 3, 31)))
        ds_TB = xr.open_mfdataset(files[time_mask], combine='nested', concat_dim='time').chunk(dict(time=-1, x=100, y=100)).where(amsr_mask)
        ds_TB['time'] = times[time_mask]
    
        # Apply melt algorithm
        melt = {}
        
        # Gapfilling and TBV mask computation
        tbh_morning = algo.gapfilling(ds_TB['TB_H_morning'].sel(iterations = 10), days=3)
        tbh_afternoon = algo.gapfilling(ds_TB['TB_H_afternoon'].sel(iterations = 10), days=3)
        tbh_mean = algo.gapfilling(ds_TB[['TB_H_afternoon', 'TB_H_morning']].sel(iterations = 10).to_array(
            dim = 'new').mean('new'), days=3)  
        mask_TBV_morning = ds_TB['TB_V_morning'].sel(iterations = 10).std(dim='time') > thres_TBV
        mask_TBV_afternoon = ds_TB['TB_V_afternoon'].sel(iterations = 10).std(dim='time') > thres_TBV
        mask_TBV_daily = ds_TB[['TB_V_afternoon', 'TB_V_morning']].sel(iterations = 10).to_array(
            dim = 'new').mean('new').std(dim='time') > thres_TBV
        
        # Compute adaptive threshold
        stats_morning = algo.compute_stats_adaptive(tbh_morning, melt_coef=3.0, threshold0=15)
        stats_afternoon = algo.compute_stats_adaptive(tbh_afternoon, melt_coef=3.0, threshold0=15)
        stats_daily = algo.compute_stats_adaptive(tbh_mean, melt_coef=3.0, threshold0=15)
        
        # Compute melt
        melt[f'{WMO_MELT_VAR}_smos_daily'] = algo.detect_melt_PF06(tbh_mean, stats_daily, melt_coef=3.0, threshold_range=(10, 25))
        melt[f'{WMO_MELT_VAR}_smos_morning'] = algo.detect_melt_PF06(tbh_morning, stats_morning, melt_coef=3.0, threshold_range=(10, 25))
        melt[f'{WMO_MELT_VAR}_smos_afternoon'] = algo.detect_melt_PF06(tbh_afternoon, stats_afternoon, melt_coef=3.0, threshold_range=(10, 25))
        melt['mask_TBV_morning'] = mask_TBV_morning
        melt['mask_TBV_afternoon'] = mask_TBV_afternoon
        melt['mask_TBV_daily'] = mask_TBV_daily
    	
    	# Export dataset
        output = xr.Dataset(melt)#.where(mask_TBV_morning, algo.missing_value)
        output.to_netcdf(path_data + "/SMOS_data/melt/%s/SMOS_melt_12km_%s_%s-%s.nc"%(name, name, y, y + 1))


def compute_contamination_mask(melt):    # Years should be a sequence containing all the m
    years = np.unique(melt.time.values.astype('datetime64[Y]').astype(int) + 1970)
    winter_melt = 0
    ndays = 0
    for y in years[:-1]:
        m = (melt > 0).sel(time = slice(f'{y}-06-15', f'{y}-08-15')).sum(['time'])
        winter_melt = winter_melt + m   
        ndays = ndays + len(melt.sel(time = slice(f'{y}-06-15', f'{y}-08-15')).time)  
    return winter_melt / ndays, years


now = time.time()
name = os.environ["NAME"]
thres_TBV = 5

smos_melt_loop('/bettik/PROJECTS/pr-snowem/zeigerp', name = name)

print('Total time ellapsed: %.2f hours'%((time.time()-now)/3600))
    
    
    


