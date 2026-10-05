# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:19:54 2023

@author: Pierre
"""

import numpy as np
import xarray as xr
import scipy.io as sio
import glob
import os
import time
import numba
#from tqdm import tqdm
from pyproj import Transformer
from scipy.optimize import curve_fit
from datetime import datetime
from multiprocessing import Pool
from multifresnel import compute_multi_alternating_layer_reflection2
from rSIR import extract_MRFs,AVE_direct,rSIR_iter
from edit_nc import edit_nc_attributes
from collections.abc import Iterable

# Correct paths and options
platform = "local"
if platform=="local":
    name = 'Antarctica'
    
elif platform=="dahu":
    name = os.environ["NAME"]
    
# Code parameters
res = 12500                 # Grid spacing in meters 
# single_incidence = True     # Whether to compute single or multi-incidence products
separate_asc = True         
max_obs_per_DGG = 20
bin_sep = 2
thres_rad_acc = 5 #2.8 
#max_incidence = 40 #55
min_pix_per_MRF = 2
max_dist_pixel = 12.5   # In km
max_iter = 11
iter_out = [0, 10] #np.arange(0,16,5)
list_incidences = [22.5, 40]    #np.arange(2.5, 40, 5).tolist() + [40] + np.arange(42.5, 65, 5).tolist()   #[12.5, 22.5, 40]
# normalization_angle = list_incidences   #52.5
# diagnostic_fit = False


def mat_to_xr(data_loc,data_TB,times,subset=None,name='Antarctica'):
    # Organize file variables
    DGG_ID,lats,lons = np.concatenate([np.repeat(data_loc[i,:3][:,np.newaxis],len(data_TB[i,0]),axis=1) for i in range(len(data_loc))],axis=1)
    
    if data_TB[0,0].shape[1] == 10:     # New "lighter" mat version with only 10 essential variables
        flags,incidence,azimuth,snapshot_ID,a,b,TB_H,TB_V,TB_H_acc,TB_V_acc = np.transpose(np.concatenate([data_TB[i,0] for i in range(len(data_TB))]))

    elif data_TB[0,0].shape[1] > 10:    # All version from Philippe, assumes the 19 variables are present in mat files
        flags,_,_,_,incidence,azimuth,_,_,snapshot_ID,a,b,TB_H,TB_V,_,_,TB_H_acc,TB_V_acc,_,_ = np.transpose(np.concatenate(
            [data_TB[i,0] for i in range(len(data_TB))]))
        
    if name=='Antarctica':
        proj = Transformer.from_crs(4326,3976,always_xy=True)
        [X,Y] = proj.transform(lons,lats)
    elif name=='Greenland':
        proj = Transformer.from_crs(4326,3413,always_xy=True)
        [X,Y] = proj.transform(lons,lats)
    else:
        raise ValueError("Invalid value for 'name': must be 'Antarctica' or 'Greenland'")
    hours = solar_time(lons,times)
    hours[hours<0] = hours[hours<0]+24
    hours[hours>=24] = hours[hours>=24]-24
    
    if subset is not None:      # Select subset
        mask = ((X>=subset[0]) & (X<=subset[1]) & (Y>=subset[2]) & (Y<=subset[3]))
    else:                       # All the observations
        mask = np.ones(len(X)).astype(bool)
    
    if len(X[mask])>0:  # Compute flags and duplicate masks for the dataset
        mask_duplicate = unique_DGG_snapshot(DGG_ID[mask],snapshot_ID[mask])
        mask_BT = (~np.isnan(TB_H[mask]) & ~np.isnan(TB_V[mask]))
        mask_flags = flag_out(flags[mask])
        mask[mask] = (mask_flags & mask_duplicate & mask_BT)
    
    ds_mat = {'DGG_ID':DGG_ID[mask],'incidence':incidence[mask],'azimuth':azimuth[mask],
              'minor_axis':b[mask],'major_axis':a[mask],'TB_H':TB_H[mask],'TB_V':TB_V[mask],
              'TB_H_acc':TB_H_acc[mask],'TB_V_acc':TB_V_acc[mask],'x':X[mask],'y':Y[mask],'hours':hours[mask]}
    # else:
    #     ds_mat = None
    #     mask = []
    return ds_mat,len(mask),np.count_nonzero(mask)
    
    
def unique_DGG_snapshot(DGG_ID,snapshot): #eliminate duplicated combinations of DGG_ID and snapshot_ID
    DGG_snap = np.asarray(["%.0f%.0f"%(DGG,snapshot[i]) for i,DGG in enumerate(DGG_ID)])
    mask = np.concatenate((np.diff(DGG_snap.astype(float)),[1]))!=0
    return mask


def extract_variables(dict_in):
    DGG_ID = np.concatenate([dict_in[i]["DGG_ID"] for i in range(len(dict_in))])
    BT_H = np.concatenate([dict_in[i]["TB_H"] for i in range(len(dict_in))])#[mask_duplicate]
    BT_V = np.concatenate([dict_in[i]["TB_V"] for i in range(len(dict_in))])#[mask_duplicate]
    X = np.concatenate([dict_in[i]["x"] for i in range(len(dict_in))])#[mask_all] #* 1000
    Y = np.concatenate([dict_in[i]["y"] for i in range(len(dict_in))])#[mask_all] #* 1000
    a = np.concatenate([dict_in[i]["major_axis"] for i in range(len(dict_in))]) * 1000 #[mask_all]
    b = np.concatenate([dict_in[i]["minor_axis"] for i in range(len(dict_in))]) * 1000 #[mask_all] 
    BT_H_acc = np.concatenate([dict_in[i]["TB_H_acc"] for i in range(len(dict_in))])#[mask_all]
    BT_V_acc = np.concatenate([dict_in[i]["TB_V_acc"] for i in range(len(dict_in))])#[mask_all]
    azimuth = np.concatenate([dict_in[i]["azimuth"] for i in range(len(dict_in))])#[mask_all]
    incidence = np.concatenate([dict_in[i]["incidence"] for i in range(len(dict_in))])#[mask_all]
    hours = np.concatenate([dict_in[i]["hours"] for i in range(len(dict_in))])#[mask_all]
    
    return np.asarray([X,Y,a,b,BT_H,BT_V,BT_H_acc,BT_V_acc,azimuth,incidence,DGG_ID,hours])
    

def separate_ascendant(datas):
    mask_none = np.asarray([d is not None for d in datas])
    print(mask_none)
    variables = extract_variables(np.asarray(datas)[mask_none])
    variables_asc = variables[:,variables[-1]<12]       # Ascendant tracks are morning tracks
    variables_desc = variables[:,variables[-1]>=12]     # Descendant tracks are afternoon tracks
    
    return variables_asc,variables_desc


def solar_time(lons,times):
    return times.hour + times.minute/60 + lons/15
    

def read_mat_day(args):
    assert isinstance(list_incidences, Iterable), "list_incidences should be iterable instance"
    assert len(list_incidences) >= 1, "list_incidences should have length >= 1"
        # Variables
    print("Into read_mat_day function")
    path, grid, subset = args
    
    if subset is not None:     
        print(subset)
        # The grid might be split in two East-West domains with overlap
        assert type(subset) == str, "subset must be a string or None"
        if subset=='West':
            filename = name + subset
            subset = [-3000000, 100000, -2500000, 2600000]
            grid = grid.where(grid.x<=50000, drop=True)
        elif subset=='East':
            filename = name + subset
            subset = [-100000, 2700000, -2500000, 2600000]
            grid = grid.where(grid.x>=-50000, drop=True)
        else:                       # All the observations
            subset = None
            filename = name
            
    else:
        filename = name
        
    x_ref = grid.x.values
    y_ref = grid.y.values
    files = glob.glob(path+'/*_%s.mat'%name)
    # Name confusion: Antarctica or Antartica
    if (len(files) == 0) & (name=='Antarctica'):
        files = glob.glob(path+'/*_%s.mat'%'Antartica')        
    
    date = path.split('/')[-3:]  #[-4:-1]
    ascendant = []
    datas = []
    lens_data = []
    print("\nReading .mat files for %s-%s-%s..."%(date[0],date[1],date[2]))
    time0 = time.time()
    for i,f in enumerate(files):
        try:
            datetime_in = datetime.strptime(f.split('/')[-1].split('_')[4],"%Y%m%dT%H%M%S")
            datetime_out = datetime.strptime(f.split('/')[-1].split('_')[5],"%Y%m%dT%H%M%S")
            fmat = sio.loadmat(f)
            data_loc = fmat['l1cdata'][0,0]   # 5 N_DGG rows, columns: DGG ID, latitude, longitude, Altitude and LandSeaMask value
            data_BT = fmat['l1cdata'][0,1]    
            ds_data,len1,len2 = mat_to_xr(data_loc,data_BT,times=datetime.fromtimestamp(sum(map(datetime.timestamp,[datetime_in,datetime_out]))/2),subset=subset,name=name)
            lens_data.append([len1,len2])
            datas.append(ds_data)
            ascendant.append(fmat['Asc'])
        except Exception as e:
            print("Invalid input: %s on file %s"%(e,f))
    lens_data = np.sum(lens_data,axis=0)
    
    print("... Done: %.2f hours ellapsed"%((time.time()-time0)/3600))
    # Select ascending and descending tracks
    data_asc, data_desc = separate_ascendant(datas)
    del datas, ds_data, data_BT, data_loc		# Remove variables to save memory
    len_asc_init, len_desc_init = [len(data_asc[0]), len(data_desc[0])]     # Count the number of morning and afternoon obs
    
    TB_H_rSIR_asc = []
    TB_V_rSIR_asc = []
    TB_H_rSIR_desc = []
    TB_V_rSIR_desc = []
    TB_H_rSIR_asc_acc = []
    TB_V_rSIR_asc_acc = []
    TB_H_rSIR_desc_acc = []
    TB_V_rSIR_desc_acc = []
    count_H_asc = []
    count_V_asc = []
    count_H_desc = []
    count_V_desc = []
    
    for incidence in list_incidences:   # Compute the MRF, perform rSIR iterations and save dataset
        print("Computing incidence %s°..."%incidence)
        min_incidence = incidence - 2.5
        max_incidence = incidence + 2.5
        time0 = time.time()
        
        # Filter out observations with radiometric accuracy greater than an absolute threshold
        variables_asc = mask_rad_acc_threshold(data_asc, thres = thres_rad_acc, max_theta = max_incidence, min_theta = min_incidence)
        variables_desc = mask_rad_acc_threshold(data_desc, thres = thres_rad_acc, max_theta = max_incidence, min_theta = min_incidence)
        # Select n_obs observations with best radiometric accuracy for each DGG
        mask_asc_H, mask_asc_V = mask_rad_acc_DGG(variables_asc, max_num = max_obs_per_DGG, max_thres = thres_rad_acc,
                                                  max_theta = max_incidence, min_theta = min_incidence)
        mask_desc_H, mask_desc_V = mask_rad_acc_DGG(variables_desc, max_num = max_obs_per_DGG, max_thres = thres_rad_acc, 
                                                    max_theta = max_incidence, min_theta = min_incidence)
        mask_asc_HV = mask_asc_H | mask_asc_V
        mask_desc_HV = mask_desc_H | mask_desc_V
        variables_asc, variables_desc = variables_asc[:, mask_asc_HV], variables_desc[:, mask_desc_HV]
        TB_H_corr_asc, TB_V_corr_asc = variables_asc[4], variables_asc[5]
        TB_H_corr_desc, TB_V_corr_desc = variables_desc[4], variables_desc[5]
        DGG_unique = np.unique(np.concatenate((variables_asc[10],variables_desc[10])))
        # DGG_vars = np.zeros((len(DGG_unique),1))
        
        # Compute MRFs
        # print("\nComputing the MRFs for all measurements...")
        MRFs_asc = extract_MRFs(x_ref, y_ref, variables_asc[0], variables_asc[1], variables_asc[2], variables_asc[3],
                                variables_asc[8], resolution = res)
        # print("... MRFs ascending tracks: %.2f hours ellapsed"%((time1-time0)/3600))
        MRFs_desc = extract_MRFs(x_ref, y_ref, variables_desc[0], variables_desc[1], variables_desc[2], variables_desc[3],
                                 variables_desc[8], resolution = res)
        # print("... MRFs descending tracks: %.2f hours ellapsed"%((time.time()-time1)/3600))
        mask_MRF_asc = np.asarray([len(MRFs_asc[i][0]) for i in range(len(MRFs_asc))]) >= min_pix_per_MRF  # Quality control for the minimal number of pixels per MRF
        mask_MRF_desc = np.asarray([len(MRFs_desc[i][0]) for i in range(len(MRFs_desc))]) >= min_pix_per_MRF  # Quality control for the minimal number of pixels per MRF
        mask_asc_H = mask_asc_H[mask_asc_HV] & (~np.isnan(TB_H_corr_asc))
        mask_asc_V = mask_asc_V[mask_asc_HV] & (~np.isnan(TB_V_corr_asc))
        mask_desc_H = mask_desc_H[mask_desc_HV] & (~np.isnan(TB_H_corr_desc))
        mask_desc_V = mask_desc_V[mask_desc_HV] & (~np.isnan(TB_V_corr_desc))
        len_H_asc_fin, len_H_desc_fin = [len(TB_H_corr_asc[mask_MRF_asc]), len(TB_H_corr_desc[mask_MRF_desc])]
        len_V_asc_fin, len_V_desc_fin = [len(TB_V_corr_asc[mask_MRF_asc]), len(TB_V_corr_desc[mask_MRF_desc])]
        
        # Compute AVE (first iteration of rSIR)
        AVE_asc = np.asarray(AVE_direct(x_ref, y_ref, variables_asc[0][mask_MRF_asc], variables_asc[1][mask_MRF_asc],
                                        MRFs_asc[mask_MRF_asc], TB_H_corr_asc[mask_MRF_asc], TB_V_corr_asc[mask_MRF_asc],
                                        variables_asc[6][mask_MRF_asc], variables_asc[7][mask_MRF_asc],
                                        mask_H_acc = mask_asc_H[mask_MRF_asc], mask_V_acc = mask_asc_V[mask_MRF_asc], dist_pix = max_dist_pixel))
        AVE_desc = np.asarray(AVE_direct(x_ref, y_ref, variables_desc[0][mask_MRF_desc], variables_desc[1][mask_MRF_desc],
                                         MRFs_desc[mask_MRF_desc], TB_H_corr_desc[mask_MRF_desc], TB_V_corr_desc[mask_MRF_desc],
                                         variables_desc[6][mask_MRF_desc], variables_desc[7][mask_MRF_desc],
                                         mask_H_acc = mask_desc_H[mask_MRF_desc], mask_V_acc = mask_desc_V[mask_MRF_desc], dist_pix = max_dist_pixel))
        # print("... AVE descending tracks: %.2f hours ellapsed"%((time.time()-time1)/3600))
                                     
        # Perform rSIR iterations
        TB_H_rSIR_asc.append(rSIR_iter(MRFs_asc[mask_MRF_asc], TB_H_corr_asc[mask_MRF_asc], AVE_asc[0], max_iter = max_iter, mask_acc= mask_asc_H[mask_MRF_asc]))
        TB_V_rSIR_asc.append(rSIR_iter(MRFs_asc[mask_MRF_asc], TB_V_corr_asc[mask_MRF_asc], AVE_asc[1], max_iter = max_iter, mask_acc = mask_asc_V[mask_MRF_asc]))
        TB_H_rSIR_desc.append(rSIR_iter(MRFs_desc[mask_MRF_desc], TB_H_corr_desc[mask_MRF_desc], AVE_desc[0], max_iter = max_iter, mask_acc = mask_desc_H[mask_MRF_desc]))
        TB_V_rSIR_desc.append(rSIR_iter(MRFs_desc[mask_MRF_desc], TB_V_corr_desc[mask_MRF_desc], AVE_desc[1], max_iter = max_iter, mask_acc = mask_desc_V[mask_MRF_desc]))
    
        # Add radiometric accuracies at first iteration (weighted average) + L1C obsrvation count
        TB_H_rSIR_asc_acc.append(AVE_asc[2])
        TB_V_rSIR_asc_acc.append(AVE_asc[3])
        TB_H_rSIR_desc_acc.append(AVE_desc[2])
        TB_V_rSIR_desc_acc.append(AVE_desc[3])
        count_H_asc.append(AVE_asc[4])
        count_V_asc.append(AVE_asc[5])
        count_H_desc.append(AVE_desc[4])
        count_V_desc.append(AVE_desc[5])
        
        print("Incidence %s°: %.2f hours ellapsed"%(incidence, (time.time() - time0) / 3600))
    
    # Reshape TB datasets with four dimensions: incidence, iterations, y and x
    TB_H_rSIR_asc = np.array(TB_H_rSIR_asc).reshape(
        len(list_incidences), max_iter + 1, len(y_ref), len(x_ref))[:, iter_out, :, :]
    TB_V_rSIR_asc = np.array(TB_V_rSIR_asc).reshape(
        len(list_incidences), max_iter + 1, len(y_ref), len(x_ref))[:, iter_out, :, :]
    TB_H_rSIR_desc = np.array(TB_H_rSIR_desc).reshape(
        len(list_incidences), max_iter + 1, len(y_ref), len(x_ref))[:, iter_out, :, :]
    TB_V_rSIR_desc = np.array(TB_V_rSIR_desc).reshape(
        len(list_incidences), max_iter + 1, len(y_ref), len(x_ref))[:, iter_out, :, :]
    
    # Idem for radiometric accuracies (weighted average) and pixel counts with dimensions: incidence, y and x
    TB_H_rSIR_asc_acc = np.array(TB_H_rSIR_asc_acc).reshape(len(list_incidences), len(y_ref), len(x_ref))
    TB_V_rSIR_asc_acc = np.array(TB_V_rSIR_asc_acc).reshape(len(list_incidences), len(y_ref), len(x_ref))
    TB_H_rSIR_desc_acc = np.array(TB_H_rSIR_desc_acc).reshape(len(list_incidences), len(y_ref), len(x_ref))
    TB_V_rSIR_desc_acc = np.array(TB_V_rSIR_desc_acc).reshape(len(list_incidences), len(y_ref), len(x_ref))
    count_H_asc = np.array(count_H_asc).reshape(len(list_incidences), len(y_ref), len(x_ref))
    count_V_asc = np.array(count_V_asc).reshape(len(list_incidences), len(y_ref), len(x_ref))
    count_H_desc = np.array(count_H_desc).reshape(len(list_incidences), len(y_ref), len(x_ref))
    count_V_desc = np.array(count_V_desc).reshape(len(list_incidences), len(y_ref), len(x_ref))
    
    # Create and export output dataset
    ds_TB = xr.Dataset({"TB_H_morning":xr.DataArray(TB_H_rSIR_asc, coords=[list_incidences, iter_out, y_ref, x_ref],
                                                    dims=["incidence", "iterations", "y", "x"]),
                        "TB_V_morning":xr.DataArray(TB_V_rSIR_asc, coords=[list_incidences, iter_out, y_ref, x_ref],
                                                    dims=["incidence", "iterations", "y", "x"]),
                        "TB_H_afternoon":xr.DataArray(TB_H_rSIR_desc, coords=[list_incidences, iter_out, y_ref, x_ref],
                                                    dims=["incidence", "iterations", "y", "x"]),
                        "TB_V_afternoon":xr.DataArray(TB_V_rSIR_desc, coords=[list_incidences, iter_out, y_ref, x_ref],
                                                    dims=["incidence", "iterations", "y", "x"]),
                        "TB_H_acc_morning":xr.DataArray(TB_H_rSIR_asc_acc, coords=[list_incidences, y_ref, x_ref],
                                                        dims=["incidence", "y", "x"]),
                        "TB_V_acc_morning":xr.DataArray(TB_V_rSIR_asc_acc, coords=[list_incidences, y_ref, x_ref],
                                                        dims=["incidence", "y", "x"]),
                        "TB_H_acc_afternoon":xr.DataArray(TB_H_rSIR_desc_acc, coords=[list_incidences, y_ref, x_ref],
                                                        dims=["incidence", "y", "x"]),
                        "TB_V_acc_afternoon":xr.DataArray(TB_V_rSIR_desc_acc, coords=[list_incidences, y_ref, x_ref],
                                                        dims=["incidence", "y", "x"]),
                        "count_H_morning":xr.DataArray(count_H_asc, coords=[list_incidences, y_ref, x_ref],
                                                        dims=["incidence", "y", "x"]),
                        "count_H_afternoon":xr.DataArray(count_H_desc, coords=[list_incidences, y_ref, x_ref],
                                                        dims=["incidence", "y", "x"]),
                        "count_V_morning":xr.DataArray(count_V_asc, coords=[list_incidences, y_ref, x_ref],
                                                        dims=["incidence", "y", "x"]),
                        "count_V_afternoon":xr.DataArray(count_V_desc, coords=[list_incidences, y_ref, x_ref],
                                                        dims=["incidence", "y", "x"])
                        })

    ds_TB = edit_nc_attributes(ds_TB, x_ref, y_ref, lens_data, len_asc_init, len_desc_init, len(DGG_unique),
                               len_H_asc_fin, len_H_desc_fin, len_V_asc_fin, len_V_desc_fin)
    
    ds_TB.to_netcdf(path + '/SMOS_TB_rSIR_%s_%s-%s-%s.nc'%(filename, date[0], date[1], date[2]))
    
    
def reconstruct_from_mat(path_data, grid, subset=None):
    dirs = os.environ["PARAM"]
    dir_days = glob.glob(dirs+'/*')
    print(dir_days)
    p = Pool(len(dir_days))      
    p.map(read_mat_day,[[d,grid,subset] for d in dir_days])


def flag_out(flags):
    flag_bin = np.asarray(['{0:016b}'.format(f) for f in flags.astype(int)])
    sun_fov = np.asarray([f[-3]=='1' for f in flag_bin])
    sun_glint_fov = np.asarray([f[-4]=='1' for f in flag_bin])
    moon_fov = np.asarray([f[-5]=='1' for f in flag_bin])
    single_snapshot = np.asarray([f[-6]=='1' for f in flag_bin])
    sun_point = np.asarray([f[8]=='1' for f in flag_bin])
    sun_glint_area = np.asarray([f[7]=='1' for f in flag_bin])
    moon_point = np.asarray([f[6]=='1' for f in flag_bin])
    af_fov = np.asarray([f[5]=='1' for f in flag_bin])
    border_fov = np.asarray([f[3]=='1' for f in flag_bin])
    sun_tails = np.asarray([f[2]=='1' for f in flag_bin])
    rfi = np.asarray([f[9]=='1' for f in flag_bin])
    rfi_tails = np.asarray([f[4]=='1' for f in flag_bin])
    rfi_amplitude = np.asarray([f[:2]=='00' for f in flag_bin])
    # Filter out flagged observations with the least conservative option
    mask_flag = (~sun_point & ~sun_glint_area & af_fov & ~border_fov & ~rfi & ~rfi_tails & rfi_amplitude)
    return mask_flag


# Extract for each DGG the n measurements with lowest radiometric accuracy (with absolute accuracy threshold)
def mask_rad_acc_DGG(variables, max_num = 20, max_thres = 2.8, max_theta = 55, min_theta = 0):
    print("\nCompute radiometric accuracy mask in each DGG")
    max_num = int(max_num)
    mask_nan_H = (~np.isnan(variables[4])) & (variables[9]<=max_theta) & (variables[9]>=min_theta) & (variables[6]<=max_thres)
    mask_nan_V = (~np.isnan(variables[5])) & (variables[9]<=max_theta) & (variables[9]>=min_theta) & (variables[7]<=max_thres)
    TB_H_acc = variables[6][mask_nan_H]
    TB_V_acc = variables[7][mask_nan_V]
    DGG_H = variables[10][mask_nan_H]
    DGG_V = variables[10][mask_nan_V]
    DGGs = np.unique(variables[10][mask_nan_H | mask_nan_V])
    maskH, maskV = np.ones(len(DGG_H)).astype(bool), np.ones(len(DGG_V)).astype(bool)
    # It is faster this way than using 2D indexing, as the combined HV array is ~1.5 times larger than the separated H and V arrays (because radiometric accuracy differs)
    for i,dg in enumerate(DGGs):
        mask_DGG_H = DGG_H==dg
        mask_DGG_V = DGG_V==dg
        # sortedH,sortedV = np.sort(TB_H_acc[mask_DGG_H]),np.sort(TB_V_acc[mask_DGG_V])
        if np.count_nonzero(mask_DGG_H)>max_num:
            maskH[mask_DGG_H & (TB_H_acc>=np.sort(TB_H_acc[mask_DGG_H])[max_num])] = False
            # maskH[mask_DGG_H & (TB_H_acc>=TB_H_acc[mask_DGG_H][np.argpartition(TB_H_acc[mask_DGG_H],max_num)[max_num]])] = False
        if np.count_nonzero(mask_DGG_V)>max_num:
            maskV[mask_DGG_V & (TB_V_acc>=np.sort(TB_V_acc[mask_DGG_V])[max_num])] = False
            # maskV[mask_DGG_V & (TB_V_acc>=TB_V_acc[mask_DGG_V][np.argpartition(TB_V_acc[mask_DGG_V],max_num)[max_num]])] = False
    output_H, output_V = np.zeros(len(variables[10])).astype(bool), np.zeros(len(variables[10])).astype(bool)
    output_H[mask_nan_H] = maskH
    output_V[mask_nan_V] = maskV
    return output_H, output_V


# Extract for each DGG the n measurements with lowest radiometric accuracy (with absolute accuracy threshold)
def mask_rad_acc_DGG_bins(variables, max_num = 20, max_thres = 2.8, max_theta = 55, min_theta = 0):
    print("\nCompute radiometric accuracy mask in each DGG")
    max_num = int(max_num)
    mask_nan_H = (~np.isnan(variables[4])) & (variables[9]<=max_theta) & (variables[9]>=min_theta) & (variables[6]<=max_thres)
    mask_nan_V = (~np.isnan(variables[5])) & (variables[9]<=max_theta) & (variables[9]>=min_theta) & (variables[7]<=max_thres)
    TB_H_acc = variables[6][mask_nan_H]
    TB_V_acc = variables[7][mask_nan_V]
    DGG_H = variables[10][mask_nan_H]
    DGG_V = variables[10][mask_nan_V]
    incidences_H = variables[9][mask_nan_H]
    incidences_V = variables[9][mask_nan_V]
    DGGs = np.unique(variables[10][mask_nan_H | mask_nan_V])
    maskH, maskV = np.ones(len(DGG_H)).astype(bool), np.ones(len(DGG_V)).astype(bool)
    # It is faster this way than using 2D indexing, as the combined HV array is ~1.5 times larger than the separated H and V arrays (because radiometric accuracy differs)
    for i,dg in enumerate(DGGs):
        mask_DGG_H = DGG_H==dg
        mask_DGG_V = DGG_V==dg
        if np.count_nonzero(mask_DGG_H)>max_num:
            bins_H = np.arange(int(np.min(incidences_H[mask_DGG_H])), np.max(incidences_H[mask_DGG_H]), bin_sep)
            mask_incidence_H = np.isin(np.arange(np.count_nonzero(mask_DGG_H)), 
                                     rad_acc_per_bin(incidences_H[mask_DGG_H], TB_H_acc[mask_DGG_H], bins_H))
            maskH[mask_DGG_H] = mask_incidence_H
        if np.count_nonzero(mask_DGG_V)>max_num:
            bins_V = np.arange(int(np.min(incidences_V[mask_DGG_V])), np.max(incidences_V[mask_DGG_V]), bin_sep)
            mask_incidence_V = np.isin(np.arange(np.count_nonzero(mask_DGG_V)), 
                                     rad_acc_per_bin(incidences_V[mask_DGG_V], TB_V_acc[mask_DGG_V], bins_V))
            maskV[mask_DGG_V] = mask_incidence_V

    output_H, output_V = np.zeros(len(variables[10])).astype(bool), np.zeros(len(variables[10])).astype(bool)
    output_H[mask_nan_H] = maskH
    output_V[mask_nan_V] = maskV
    return output_H, output_V


@numba.jit(nopython=True, cache=True)
def rad_acc_per_bin(incidence, rad_acc, bins):
    output = []
    for b in bins:
        tmp = rad_acc[(incidence >= b) & (incidence < b + bin_sep)]
        if len(tmp) > 0:
            output.append(np.argwhere(rad_acc == np.nanmin(tmp))[0][0])
    return np.asarray(output).ravel()
        

def mask_rad_acc_threshold(variables, thres = 2.8, max_theta = 70, min_theta = 0, max_tb = 280, min_tb = [70, 110]):
    mask_nan_H = ((~np.isnan(variables[4])) & (variables[9]<=max_theta) & (variables[9]>=min_theta) & 
                  (variables[6]<=thres) & (variables[4] <= max_tb) & (variables[4] >= min_tb[0]))
    mask_nan_V = ((~np.isnan(variables[5])) & (variables[9]<=max_theta) & (variables[9]>=min_theta) & 
                  (variables[7]<=thres) & (variables[5] <= max_tb) & (variables[5] >= min_tb[1]))
    return variables[:, mask_nan_H & mask_nan_V]   # variables[:, mask_nan]


def set_DGG_variables(variables_asc, variables_desc, popts_asc, popts_desc, residu_asc, residu_desc):
    DGG_unique = np.unique(np.concatenate((variables_asc[10],variables_desc[10])))
    DGG_X = np.asarray([np.concatenate((variables_asc[0][variables_asc[10] == dg], 
                        variables_desc[0][variables_desc[10] == dg]))[0] for dg in DGG_unique]).ravel()
    DGG_Y = np.asarray([np.concatenate((variables_asc[1][variables_asc[10] == dg], 
                        variables_desc[1][variables_desc[10] == dg]))[0] for dg in DGG_unique]).ravel()
    mask_asc = np.isin(DGG_unique, np.unique(variables_asc[10]))
    mask_desc = np.isin(DGG_unique, np.unique(variables_desc[10]))
    DGG_vars = np.zeros((len(DGG_unique), 14)) * np.nan
    DGG_vars[:, 0] = DGG_X
    DGG_vars[:, 1] = DGG_Y
    DGG_vars[mask_asc, 2:7] = popts_asc
    DGG_vars[mask_asc, 7] = residu_asc
    DGG_vars[mask_desc, 8:13] = popts_desc
    DGG_vars[mask_desc, 13] = residu_desc
    
    return DGG_vars, DGG_unique


        
        






