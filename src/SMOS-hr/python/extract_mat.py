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

# Code parameters
name = os.environ["NAME"]
#name = "Greenland"
res = 12500                 # Grid spacing in meters 
single_incidence = True     # Whether to compute single or multi-incidence products
separate_asc = True         
max_obs_per_DGG = 20
bin_sep = 2
thres_rad_acc = 5 #2.8 #5
#max_incidence = 40 #55
min_pix_per_MRF = 2
max_dist_pixel = 12.5   # In km
max_iter = 15
iter_out = np.arange(0,16,5)
normalization_angle = 40 #52.5
diagnostic_fit = False


def mat_to_xr(data_loc,data_TB,times,subset=None,name='Antartica'):
    # Organize file variables
    if data_TB[0,0].shape[1]>11:
        DGG_ID,lats,lons = np.concatenate([np.repeat(data_loc[i,:3][:,np.newaxis],len(data_TB[i,0]),axis=1) for i in range(len(data_loc))],axis=1)
        flags,_,_,_,incidence,azimuth,_,_,snapshot_ID,a,b,TB_H,TB_V,_,_,TB_H_acc,TB_V_acc,_,_ = np.transpose(np.concatenate(
            [data_TB[i,0] for i in range(len(data_TB))]))
        if name=='Antartica':
            proj = Transformer.from_crs(4326,3976,always_xy=True)
            [X,Y] = proj.transform(lons,lats)
        elif name=='Greenland':
            proj = Transformer.from_crs(4326,3413,always_xy=True)
            [X,Y] = proj.transform(lons,lats)
        else:
            raise ValueError("Invalid value for 'name': must be 'Antartica' or 'Greenland'")
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
            # fig_article_SMOSHR.fig1_smos_data(snapshot_ID, incidence, X, Y, a, b, TB_H_acc, mask, indice = indice) 	# Fig 1 from RSE paper
        
        ds_mat = {'DGG_ID':DGG_ID[mask],'incidence':incidence[mask],'azimuth':azimuth[mask],
                  'minor_axis':b[mask],'major_axis':a[mask],'TB_H':TB_H[mask],'TB_V':TB_V[mask],
                  'TB_H_acc':TB_H_acc[mask],'TB_V_acc':TB_V_acc[mask],'x':X[mask],'y':Y[mask],'hours':hours[mask]}
    else:
        ds_mat = None
        mask = []
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
    # Variables
    print("Into read_mat_day function")
    path,grid,subset = args
    x_ref = grid.x.values
    y_ref = grid.y.values
    files = glob.glob(path+'/*_%s.mat'%name)
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
    if single_incidence:
        min_incidence = normalization_angle - 2.5
        max_incidence = normalization_angle + 2.5
    else:
        min_incidence = 0  
        max_incidence = 40  
    print("... Done: %.2f hours ellapsed"%((time.time()-time0)/3600))

    # Compute the MRF, perform rSIR iterations and save dataset
    if separate_asc:    # Separate computation for ascending and descending passes
        variables_asc,variables_desc = separate_ascendant(datas)
        del datas,ds_data,data_BT,data_loc		# Remove variables to save memory
        len_asc_init,len_desc_init = [len(variables_asc[0]),len(variables_desc[0])]     # Count the number of morning and afternoon obs
        
        print("\nComputing radiometric accuracy mask in each DGG...")
        time0 = time.time()
        if single_incidence:
            # Filter out observations with radiometric accuracy greater than an absolute threshold
            variables_asc = mask_rad_acc_threshold(variables_asc, thres = thres_rad_acc, max_theta = max_incidence, min_theta = min_incidence)
            variables_desc = mask_rad_acc_threshold(variables_desc, thres = thres_rad_acc, max_theta = max_incidence, min_theta = min_incidence)
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
            print("... Rad acc: %.2f hours ellapsed for morning and afternoon tracks"%((time.time()-time0)/3600))
            DGG_unique = np.unique(np.concatenate((variables_asc[10],variables_desc[10])))
            DGG_vars = np.zeros((len(DGG_unique),1))
            
        else:    
            # Filter out observations with radiometric accuracy greater than an absolute threshold
            if diagnostic_fit:
                variables_40_asc = mask_rad_acc_threshold(variables_asc, thres = thres_rad_acc, max_theta = 42.5, min_theta = 37.5)
                variables_40_desc = mask_rad_acc_threshold(variables_desc, thres = thres_rad_acc, max_theta = 42.5, min_theta = 37.5)
            variables_asc = mask_rad_acc_threshold(variables_asc, thres = thres_rad_acc, max_theta = 70, min_theta = 0)
            variables_desc = mask_rad_acc_threshold(variables_desc, thres = thres_rad_acc, max_theta = 70, min_theta = 0)
            # Select n_obs observations with best radiometric accuracy in each DGG, per incidence bin
            mask_asc_H, mask_asc_V = mask_rad_acc_DGG_bins(variables_asc, max_num = max_obs_per_DGG, max_thres = thres_rad_acc,
                                                           max_theta = 60, min_theta = min_incidence)
            mask_desc_H, mask_desc_V = mask_rad_acc_DGG_bins(variables_desc, max_num = max_obs_per_DGG, max_thres = thres_rad_acc,
                                                             max_theta = 60, min_theta = min_incidence)
            mask_asc_HV = mask_asc_H | mask_asc_V
            mask_desc_HV = mask_desc_H | mask_desc_V
            variables_asc, variables_desc = variables_asc[:, mask_asc_HV], variables_desc[:, mask_desc_HV]
            # fig_article_SMOSHR.fig3_incidence(variables_asc, variables_desc, random_choice = 7144462) 		# Fig 3 from RSE paper
            print("... Rad acc: %.2f hours ellapsed for morning and afternoon tracks"%((time.time()-time0)/3600))
            
            # Fresnel corrections for the incidence in the selected dataset
            print("\nComputing the Fresnel relationship between incidence and TB for each DGG...")
            time0 = time.time()
            TB_H_corr_asc, TB_V_corr_asc, popts_asc, residu_asc = fit_multilayer_reg(variables_asc, angle_norm = normalization_angle)     
            TB_H_corr_desc, TB_V_corr_desc, popts_desc, residu_desc = fit_multilayer_reg(variables_desc, angle_norm = normalization_angle) 
            # Export optimal parameters and residue for each DGG
            DGG_vars, DGG_unique = set_DGG_variables(variables_asc, variables_desc, popts_asc, popts_desc, residu_asc, residu_desc)
            
            if diagnostic_fit:
                export_diagnostic_fit(variables_40_asc, variables_40_desc, variables_asc, variables_desc, 
                                      TB_H_corr_asc, TB_H_corr_desc, TB_V_corr_asc, TB_V_corr_desc, DGG_vars, 
                                      filename = path + '/Diag_bias_fit_%s-%s-%s_%s.nc'%(date[0], date[1], date[2], name))
                return []
            
            # Remove high-incidence observations with low native spatial resolution
            mask_incidence_asc = variables_asc[9] <= max_incidence
            mask_incidence_desc = variables_desc[9] <= max_incidence
            TB_H_corr_asc = TB_H_corr_asc[mask_incidence_asc]
            TB_V_corr_asc = TB_V_corr_asc[mask_incidence_asc]
            TB_H_corr_desc = TB_H_corr_desc[mask_incidence_desc]
            TB_V_corr_desc = TB_V_corr_desc[mask_incidence_desc]
            mask_asc_HV[mask_asc_HV] = mask_incidence_asc
            mask_desc_HV[mask_desc_HV] = mask_incidence_desc
            variables_asc = variables_asc[:, mask_incidence_asc]
            variables_desc = variables_desc[:, mask_incidence_desc]
            print("... Fresnel: %.2f hours ellapsed for morning and afternoon tracks"%((time.time()-time0)/3600))
    
        # Compute MRFs
        print("\nComputing the MRFs for all measurements...")
        time0 = time.time()
        MRFs_asc = extract_MRFs(x_ref, y_ref, variables_asc[0], variables_asc[1], variables_asc[2], variables_asc[3],
                                variables_asc[8], resolution = res)
        time1 = time.time()
        print("... MRFs ascending tracks: %.2f hours ellapsed"%((time1-time0)/3600))
        MRFs_desc = extract_MRFs(x_ref, y_ref, variables_desc[0], variables_desc[1], variables_desc[2], variables_desc[3],
                                 variables_desc[8], resolution = res)
        print("... MRFs descending tracks: %.2f hours ellapsed"%((time.time()-time1)/3600))
        mask_MRF_asc = np.asarray([len(MRFs_asc[i][0]) for i in range(len(MRFs_asc))]) >= min_pix_per_MRF  # Quality control for the minimal number of pixels per MRF
        mask_MRF_desc = np.asarray([len(MRFs_desc[i][0]) for i in range(len(MRFs_desc))]) >= min_pix_per_MRF  # Quality control for the minimal number of pixels per MRF
        mask_asc_H = mask_asc_H[mask_asc_HV] & (~np.isnan(TB_H_corr_asc))
        mask_asc_V = mask_asc_V[mask_asc_HV] & (~np.isnan(TB_V_corr_asc))
        mask_desc_H = mask_desc_H[mask_desc_HV] & (~np.isnan(TB_H_corr_desc))
        mask_desc_V = mask_desc_V[mask_desc_HV] & (~np.isnan(TB_V_corr_desc))
        len_H_asc_fin, len_H_desc_fin = [len(TB_H_corr_asc[mask_MRF_asc]), len(TB_H_corr_desc[mask_MRF_desc])]
        len_V_asc_fin, len_V_desc_fin = [len(TB_V_corr_asc[mask_MRF_asc]), len(TB_V_corr_desc[mask_MRF_desc])]
        
        # Compute AVE (first iteration of rSIR)
        print("\nComputing rSIR first iteration (AVE)...")
        time0 = time.time()
        AVE_asc = np.asarray(AVE_direct(x_ref, y_ref, variables_asc[0][mask_MRF_asc], variables_asc[1][mask_MRF_asc],
                                        MRFs_asc[mask_MRF_asc], TB_H_corr_asc[mask_MRF_asc], TB_V_corr_asc[mask_MRF_asc],
                                        variables_asc[6][mask_MRF_asc], variables_asc[7][mask_MRF_asc],
                                        mask_H_acc = mask_asc_H[mask_MRF_asc], mask_V_acc = mask_asc_V[mask_MRF_asc], dist_pix = max_dist_pixel))
        time1 = time.time()
        print("... AVE ascending tracks: %.2f hours ellapsed"%((time1-time0)/3600))
        AVE_desc = np.asarray(AVE_direct(x_ref, y_ref, variables_desc[0][mask_MRF_desc], variables_desc[1][mask_MRF_desc],
                                         MRFs_desc[mask_MRF_desc], TB_H_corr_desc[mask_MRF_desc], TB_V_corr_desc[mask_MRF_desc],
                                         variables_desc[6][mask_MRF_desc], variables_desc[7][mask_MRF_desc],
                                         mask_H_acc = mask_desc_H[mask_MRF_desc], mask_V_acc = mask_desc_V[mask_MRF_desc], dist_pix = max_dist_pixel))
        print("... AVE descending tracks: %.2f hours ellapsed"%((time.time()-time1)/3600))
                                     
        # Perform rSIR iterations
        TB_H_rSIR_asc = rSIR_iter(MRFs_asc[mask_MRF_asc], TB_H_corr_asc[mask_MRF_asc], AVE_asc[0], max_iter = max_iter, mask_acc= mask_asc_H[mask_MRF_asc])
        TB_V_rSIR_asc = rSIR_iter(MRFs_asc[mask_MRF_asc], TB_V_corr_asc[mask_MRF_asc], AVE_asc[1], max_iter = max_iter, mask_acc = mask_asc_V[mask_MRF_asc])
        TB_H_rSIR_desc = rSIR_iter(MRFs_desc[mask_MRF_desc], TB_H_corr_desc[mask_MRF_desc], AVE_desc[0], max_iter = max_iter, mask_acc = mask_desc_H[mask_MRF_desc])
        TB_V_rSIR_desc = rSIR_iter(MRFs_desc[mask_MRF_desc], TB_V_corr_desc[mask_MRF_desc], AVE_desc[1], max_iter = max_iter, mask_acc = mask_desc_V[mask_MRF_desc])
        
        # Create and export dataset
        ds_TB = xr.Dataset({"TB_H_morning":xr.DataArray(TB_H_rSIR_asc.reshape(TB_H_rSIR_asc.shape[0],len(y_ref),len(x_ref))[iter_out],
                                                    coords=[iter_out,y_ref,x_ref],dims=["iterations","y","x"]),
                            "TB_V_morning":xr.DataArray(TB_V_rSIR_asc.reshape(TB_V_rSIR_asc.shape[0],len(y_ref),len(x_ref))[iter_out],
                                                    coords=[iter_out,y_ref,x_ref],dims=["iterations","y","x"]),
                            "TB_H_afternoon":xr.DataArray(TB_H_rSIR_desc.reshape(TB_H_rSIR_desc.shape[0],len(y_ref),len(x_ref))[iter_out],
                                                    coords=[iter_out,y_ref,x_ref],dims=["iterations","y","x"]),
                            "TB_V_afternoon":xr.DataArray(TB_V_rSIR_desc.reshape(TB_V_rSIR_desc.shape[0],len(y_ref),len(x_ref))[iter_out],
                                                    coords=[iter_out,y_ref,x_ref],dims=["iterations","y","x"]),
                            "TB_H_acc_morning":xr.DataArray(AVE_asc[2].reshape(len(y_ref),len(x_ref)),coords=[y_ref,x_ref],dims=["y","x"]),
                            "TB_V_acc_morning":xr.DataArray(AVE_asc[3].reshape(len(y_ref),len(x_ref)),coords=[y_ref,x_ref],dims=["y","x"]),
                            "TB_H_acc_afternoon":xr.DataArray(AVE_desc[2].reshape(len(y_ref),len(x_ref)),coords=[y_ref,x_ref],dims=["y","x"]),
                            "TB_V_acc_afternoon":xr.DataArray(AVE_desc[3].reshape(len(y_ref),len(x_ref)),coords=[y_ref,x_ref],dims=["y","x"]),
                            "count_H_morning":xr.DataArray(AVE_asc[4].reshape(len(y_ref),len(x_ref)),coords=[y_ref,x_ref],dims=["y","x"]),
                            "count_H_afternoon":xr.DataArray(AVE_desc[4].reshape(len(y_ref),len(x_ref)),coords=[y_ref,x_ref],dims=["y","x"]),
                            "count_V_morning":xr.DataArray(AVE_asc[5].reshape(len(y_ref),len(x_ref)),coords=[y_ref,x_ref],dims=["y","x"]),
                            "count_V_afternoon":xr.DataArray(AVE_desc[5].reshape(len(y_ref),len(x_ref)),coords=[y_ref,x_ref],dims=["y","x"])})
        if single_incidence==False:
            ds_TB["multilayer_fit_params"] = xr.DataArray(DGG_vars, coords = [DGG_unique, np.arange(len(DGG_vars[0]))], dims = ["DGG_ID", "variables"])
        ds_TB = edit_nc_attributes(ds_TB, x_ref, y_ref, lens_data, len_asc_init, len_desc_init, len(DGG_unique),
                                   len_H_asc_fin, len_H_desc_fin, len_V_asc_fin, len_V_desc_fin, single_incidence)
        
        if single_incidence:
            ds_TB.to_netcdf(path + '/TB_rSIR_%s_%s-%s-%s_single_incidence_%s.nc'%(name, date[0], date[1], date[2], normalization_angle))
        else:
            ds_TB.to_netcdf(path + '/TB_rSIR_%s_%s-%s-%s_low-incidence-%s.nc'%(name, date[0], date[1], date[2], max_incidence))
    
    
def reconstruct_from_mat(path_data,grid,subset=None):
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


def fit_multilayer_reg(variables, DGGs = None, angle_norm = 52.5, max_num = max_obs_per_DGG):   # Fit dielectric constants, T, alpha and beta from ascending & descending H and V measurements
    incidences = variables[9]
    TB_H = variables[4]
    TB_V = variables[5]
    DGG = variables[10]
    mus = np.cos(incidences * np.pi / 180)
    single_incidence = (incidences >= angle_norm - 2.5) & (incidences <= angle_norm + 2.5)
    TB_H_corr = np.zeros(len(TB_H))*np.nan
    TB_V_corr = np.zeros(len(TB_V))*np.nan
    TB_H_corr[single_incidence] = TB_H[single_incidence]
    TB_V_corr[single_incidence] = TB_V[single_incidence]
    
    print("\nComputing the relationship between incidence and TB for each DGG...")
    if DGGs is None:
        DGGs = np.unique(DGG)
    coefs = []
    residu = []
    first_guess_eps2 = 1.78
    first_guess_eps3 = 1.92
    first_guess_alpha = 0.25
    first_guess_beta = 0.97
    first_guess_T = 225
    for i,dg in enumerate(DGGs):
        mask_DGG = DGG==dg
        if ((len(incidences[mask_DGG]) >= 10) & (np.nanmax(incidences[mask_DGG]) > angle_norm - 2.5)):
            try:
                popts,pcov = curve_fit(cost_function_reg, np.hstack([mus[mask_DGG], TB_V[mask_DGG], TB_H[mask_DGG]]), 
                                        np.zeros(len(mus[mask_DGG]) * 4), 
                                        p0 = [first_guess_eps2, first_guess_eps3, first_guess_T, first_guess_alpha, first_guess_beta], 
                                        bounds = [(1, 1, 0, -2, 0.8), (5, 5, 300, 2, 1.2)], xtol = 1e-5, ftol = 1e-3)
                eV, eH = fresnel_multilayer(mus[mask_DGG], popts[0], popts[1])
                eV0, eH0 = fresnel_multilayer(np.cos(angle_norm * np.pi / 180), popts[0], popts[1])
                coefs.append(popts)
                # Compute corrected TB H and V with the coefficients
                TbH = popts[2] * eH * (popts[4] + (1 - mus[mask_DGG]) * popts[3])
                TbH0 = popts[2] * eH0 * (popts[4] + (1 - np.cos(angle_norm * np.pi / 180)) * popts[3])
                TbV = popts[2] * eV
                TbV0 = popts[2] * eV0
                TB_H_corr[mask_DGG] = TB_H[mask_DGG] * (TbH0 / TbH)
                TB_V_corr[mask_DGG] = TB_V[mask_DGG] * (TbV0 / TbV)
                residu.append(np.linalg.norm(TbH - TB_H[mask_DGG]) / np.sqrt(len(TB_H[mask_DGG])))
            except RuntimeError:
                coefs.append([np.nan,np.nan,np.nan,np.nan,np.nan])
                residu.append(np.nan)
                
        else:
            coefs.append([np.nan,np.nan,np.nan,np.nan,np.nan])
            residu.append(np.nan)

    return TB_H_corr, TB_V_corr, np.asarray(coefs), np.asarray(residu)
    
           
def cost_function_reg(incidences_obs, eps2, eps3, T, alpha, beta):
    thickness = 0.3
    kd = thickness * (2 * np.pi * 1.4e9) / 3e8   # kd = thickness * k with k = 2 PI freq / C_SPEED
    npair = 128
    eps0 = 1
    mu = incidences_obs[: int(len(incidences_obs) / 3)]
    TB_V_obs = incidences_obs[int(len(incidences_obs) / 3) : int(2 * len(incidences_obs) / 3)]
    TB_H_obs = incidences_obs[int(2 * len(incidences_obs) / 3) :]
    Rv, Rh = compute_multi_alternating_layer_reflection2(eps0, eps2, eps3, npair, kd, mu)
    TB_V = T * (1 - Rv)
    TB_H = T * (1 - Rh) * (beta + (1 - mu) * alpha)
    
    # Regularization for beta and alpha to be close to 1 and 0, respectively
    reg_weight = 2  #np.sqrt(len(TB_V))         # 2 is the ~ noise level of SMOS in K
    return np.hstack([TB_V - TB_V_obs, TB_H - TB_H_obs, np.repeat(reg_weight * (beta - 1), len(TB_V)), 
                      np.repeat(reg_weight * (alpha - 0), len(TB_V))])

 
def fresnel_multilayer(mu, eps2, eps3):       # , alpha, beta
    thickness = 0.3
    kd = thickness * (2 * np.pi * 1.4e9) / 3e8   # kd = thickness * k with k = 2 PI freq / C_SPEED
    npair = 128
    eps0 = 1
    Rv, Rh = compute_multi_alternating_layer_reflection2(eps0, eps2, eps3, npair, kd, mu)
    return (1 - Rv), (1 - Rh) #* (1 + linear_termH)

    
def fit_fresnel(variables,DGGs=None,angle_norm=52.5):   # Fit dielectric constant and T from ascending & descending H and V measurements
    incidences = variables[9]
    TB_H = variables[4]
    TB_V = variables[5]
    DGG = variables[10]
    DGGs = np.unique(DGG)
    TB_H_corr = np.zeros(len(TB_H))*np.nan
    TB_V_corr = np.zeros(len(TB_V))*np.nan
    for i,dg in enumerate(DGGs):
        mask_DGG = DGG==dg
        if len(incidences[mask_DGG])>2:
            # Model relationship between incidence and TB with the Fresnel equations
            eps_opt,pcov = curve_fit(fresnel_HV,incidences[mask_DGG],TB_H[mask_DGG] / TB_V[mask_DGG],p0=[1.3])  #,sigma=TB_H_acc[mask_DGG]+TB_V_acc[mask_DGG]
            TB_H_corr[mask_DGG] = TB_H[mask_DGG] * fresnel_H(angle_norm,eps_opt,250) / fresnel_H(incidences[mask_DGG],eps_opt,250)
            TB_V_corr[mask_DGG] = TB_V[mask_DGG] * fresnel_V(angle_norm,eps_opt,250) / fresnel_V(incidences[mask_DGG],eps_opt,250)
            
    return TB_H_corr,TB_V_corr

    
def fresnel_H(theta,epsilon,T):     # Prise en compte des indices de réfraction ?
    return T * (1 - ((np.cos(theta*np.pi/180) - np.sqrt(epsilon**2 - np.sin(theta*np.pi/180)**2))/
                     (np.cos(theta*np.pi/180) + np.sqrt(epsilon**2 - np.sin(theta*np.pi/180)**2)))**2)


def fresnel_V(theta,epsilon,T):
    return T * (1 - ((epsilon**2 * np.cos(theta*np.pi/180) - np.sqrt(epsilon**2 - np.sin(theta*np.pi/180)**2))/
                     (epsilon**2 * np.cos(theta*np.pi/180) + np.sqrt(epsilon**2 - np.sin(theta*np.pi/180)**2)))**2)


def fresnel_HV(theta,epsilon):
    eh = 1 - ((np.cos(theta*np.pi/180) - np.sqrt(epsilon**2 - np.sin(theta*np.pi/180)**2))/
              (np.cos(theta*np.pi/180) + np.sqrt(epsilon**2 - np.sin(theta*np.pi/180)**2)))**2
    ev = 1 - ((epsilon**2 * np.cos(theta*np.pi/180) - np.sqrt(epsilon**2 - np.sin(theta*np.pi/180)**2))/
              (epsilon**2 * np.cos(theta*np.pi/180) + np.sqrt(epsilon**2 - np.sin(theta*np.pi/180)**2)))**2
    return eh/ev


def fresnel_T(theta_eps,T):
    theta,epsilon = theta_eps
    return fresnel_H(theta,epsilon,T)+fresnel_V(theta,epsilon,T)


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


def export_diagnostic_fit(variables_40_asc, variables_40_desc, variables_asc, variables_desc, 
                          TB_H_corr_asc, TB_H_corr_desc, TB_V_corr_asc, TB_V_corr_desc, DGG_vars, filename):
    DGG_40 = np.unique(np.concatenate((variables_40_asc[10], variables_40_desc[10])))
    DGG_fit = np.unique(np.concatenate((variables_asc[10], variables_desc[10])))
    
    bias = [[np.nanmean(TB_H_corr_asc[variables_asc[10] == d]) - np.nanmean(variables_40_asc[4][variables_40_asc[10] == d]),
             np.nanmean(TB_H_corr_desc[variables_desc[10] == d]) - np.nanmean(variables_40_desc[4][variables_40_desc[10] == d]),
             np.nanmean(TB_V_corr_asc[variables_asc[10] == d]) - np.nanmean(variables_40_asc[5][variables_40_asc[10] == d]),
             np.nanmean(TB_V_corr_desc[variables_desc[10] == d]) - np.nanmean(variables_40_desc[5][variables_40_desc[10] == d])]
            for d in DGG_fit]
    
    ds_out = xr.Dataset({
        'DGG_vars':xr.DataArray(DGG_vars, coords=(DGG_fit, np.arange(DGG_vars.shape[1])), dims = ('DGG_ID', 'variables')),
        'DGG_bias':xr.DataArray(np.asarray(bias), coords=(DGG_fit, np.arange(4)), dims = ('DGG_ID', 'H_V_asc_desc')),
        })
    ds_out.attrs['DGG_vars_dim'] = '[0: X, 1: Y, 2-6: popts_asc, 7: residue_asc, 8-12: popts_desc, 13: residue_desc'
    ds_out.attrs['H_V_asc_desc'] = '[0: TB_H_asc, 1: TB_H_desc, 2: TB_V_asc, 3: TB_V_desc]'
    ds_out.to_netcdf(filename)

        
        






