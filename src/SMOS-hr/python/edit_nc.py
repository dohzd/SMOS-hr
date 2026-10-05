# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 10:25:07 2023

Edit attributes for distributed netCDF files

@author: Pierre
"""

import xarray as xr
import numpy as np

########################################
# Edit attributes for SMOS enhanced TB #
########################################


def edit_nc_attributes(ds_TB,x_ref,y_ref,lens_data,len_asc_init,len_desc_init,len_DGG,
                       len_H_asc_fin,len_H_desc_fin,len_V_asc_fin,len_V_desc_fin,single_incidence=True):      # Write attributes to the ds_TB dataset and contained dataArrays

    # X array
    ds_TB.x.attrs['long_name'] = 'X_coordinate_NSIDC_Polar_Stereographic_South'
    ds_TB.x.attrs['standard_name'] = 'X'
    ds_TB.x.attrs['units'] = 'meters'
    ds_TB.x.attrs['valid_min'] = '-3289335.29'
    ds_TB.x.attrs['valid_max'] = '3289335.29'
    ds_TB.x.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(x_ref),np.nanmax(x_ref))
    # Y array
    ds_TB.y.attrs['long_name'] = 'Y_coordinate_NSIDC_Polar_Stereographic_South'
    ds_TB.y.attrs['standard_name'] = 'Y'
    ds_TB.y.attrs['units'] = 'meters'
    ds_TB.y.attrs['valid_min'] = '-3323160.27'
    ds_TB.y.attrs['valid_max'] = '3323160.27'
    ds_TB.y.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(y_ref),np.nanmax(y_ref))
    # Iteration array
    ds_TB.iterations.attrs['long_name'] = 'rSIR_iterations'
    ds_TB.iterations.attrs['standard_name'] = 'Iterations'
    ds_TB.iterations.attrs['units'] = 'without_units'
    ds_TB.iterations.attrs['valid_min'] = '0'
    ds_TB.iterations.attrs['valid_max'] = '9999'
    ds_TB.iterations.attrs['actual_range'] = '[0  15]'
    # Incidence array
    ds_TB.incidence.attrs['long_name'] = 'Mean_bin_incidence_angle'
    ds_TB.incidence.attrs['standard_name'] = 'Incidence_angle'
    ds_TB.incidence.attrs['units'] = 'degree'
    ds_TB.incidence.attrs['valid_min'] = '0'
    ds_TB.incidence.attrs['valid_max'] = '90'
    ds_TB.incidence.attrs['actual_range'] = '[2.5, 62.5]'
    # TB H asc array
    ds_TB.TB_H_morning.attrs['long_name'] = 'Reconstructed_TB_H_polarization_morning_tracks'
    ds_TB.TB_H_morning.attrs['standard_name'] = 'TB_H_morning'
    ds_TB.TB_H_morning.attrs['units'] = 'Kelvins'
    ds_TB.TB_H_morning.attrs['valid_min'] = '0'
    ds_TB.TB_H_morning.attrs['valid_max'] = '300'
    ds_TB.TB_H_morning.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.TB_H_morning),np.nanmax(ds_TB.TB_H_morning))
    # TB H desc array
    ds_TB.TB_H_afternoon.attrs['long_name'] = 'Reconstructed_TB_H_polarization_afternoon_tracks'
    ds_TB.TB_H_afternoon.attrs['standard_name'] = 'TB_H_afternoon'
    ds_TB.TB_H_afternoon.attrs['units'] = 'Kelvins'
    ds_TB.TB_H_afternoon.attrs['valid_min'] = '0'
    ds_TB.TB_H_afternoon.attrs['valid_max'] = '300'
    ds_TB.TB_H_afternoon.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.TB_H_afternoon),np.nanmax(ds_TB.TB_H_afternoon))
    # TB V asc array
    ds_TB.TB_V_morning.attrs['long_name'] = 'Reconstructed_TB_V_polarization_morning_tracks'
    ds_TB.TB_V_morning.attrs['standard_name'] = 'TB_V_morning'
    ds_TB.TB_V_morning.attrs['units'] = 'Kelvins'
    ds_TB.TB_V_morning.attrs['valid_min'] = '0'
    ds_TB.TB_V_morning.attrs['valid_max'] = '300'
    ds_TB.TB_V_morning.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.TB_V_morning),np.nanmax(ds_TB.TB_V_morning))
    # TB V desc array
    ds_TB.TB_V_afternoon.attrs['long_name'] = 'Reconstructed_TB_V_polarization_afternoon_tracks'
    ds_TB.TB_V_afternoon.attrs['standard_name'] = 'TB_V_afternoon'
    ds_TB.TB_V_afternoon.attrs['units'] = 'Kelvins'
    ds_TB.TB_V_afternoon.attrs['valid_min'] = '0'
    ds_TB.TB_V_afternoon.attrs['valid_max'] = '300'
    ds_TB.TB_V_afternoon.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.TB_V_afternoon),np.nanmax(ds_TB.TB_V_afternoon))
    # TB H acc asc array
    ds_TB.TB_H_acc_morning.attrs['long_name'] = 'Reconstructed_TB_H_polarization_accuracy_morning_tracks'
    ds_TB.TB_H_acc_morning.attrs['standard_name'] = 'TB_H_morning_accuracy'
    ds_TB.TB_H_acc_morning.attrs['units'] = 'Kelvins'
    ds_TB.TB_H_acc_morning.attrs['valid_min'] = '0'
    ds_TB.TB_H_acc_morning.attrs['valid_max'] = '3'
    ds_TB.TB_H_acc_morning.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.TB_H_acc_morning),np.nanmax(ds_TB.TB_H_acc_morning))
    # TB H acc desc array
    ds_TB.TB_H_acc_afternoon.attrs['long_name'] = 'Reconstructed_TB_H_polarization_accuracy_afternoon_tracks'
    ds_TB.TB_H_acc_afternoon.attrs['standard_name'] = 'TB_H_afternoon_accuracy'
    ds_TB.TB_H_acc_afternoon.attrs['units'] = 'Kelvins'
    ds_TB.TB_H_acc_afternoon.attrs['valid_min'] = '0'
    ds_TB.TB_H_acc_afternoon.attrs['valid_max'] = '3'
    ds_TB.TB_H_acc_afternoon.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.TB_H_acc_afternoon),np.nanmax(ds_TB.TB_H_acc_afternoon))
    # TB V acc asc array
    ds_TB.TB_V_acc_morning.attrs['long_name'] = 'Reconstructed_TB_V_polarization_accuracy_morning_tracks'
    ds_TB.TB_V_acc_morning.attrs['standard_name'] = 'TB_V_morning_accuracy'
    ds_TB.TB_V_acc_morning.attrs['units'] = 'Kelvins'
    ds_TB.TB_V_acc_morning.attrs['valid_min'] = '0'
    ds_TB.TB_V_acc_morning.attrs['valid_max'] = '3'
    ds_TB.TB_V_acc_morning.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.TB_V_acc_morning),np.nanmax(ds_TB.TB_V_acc_morning))
    # TB V acc desc array
    ds_TB.TB_V_acc_afternoon.attrs['long_name'] = 'Reconstructed_TB_V_polarization_accuracy_afternoon_tracks'
    ds_TB.TB_V_acc_afternoon.attrs['standard_name'] = 'TB_V_afternoon_accuracy'
    ds_TB.TB_V_acc_afternoon.attrs['units'] = 'Kelvins'
    ds_TB.TB_V_acc_afternoon.attrs['valid_min'] = '0'
    ds_TB.TB_V_acc_afternoon.attrs['valid_max'] = '3'
    ds_TB.TB_V_acc_afternoon.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.TB_V_acc_afternoon),np.nanmax(ds_TB.TB_V_acc_afternoon))
    # Count H asc array
    ds_TB.count_H_morning.attrs['long_name'] = 'Number of SMOS observations per pixel used in TB H polarization morning reconstruction'
    ds_TB.count_H_morning.attrs['standard_name'] = 'Number of H-pol morning obs per pixel'
    ds_TB.count_H_morning.attrs['units'] = 'Kelvins'
    ds_TB.count_H_morning.attrs['valid_min'] = '0'
    ds_TB.count_H_morning.attrs['valid_max'] = '9999'
    ds_TB.count_H_morning.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.count_H_morning),np.nanmax(ds_TB.count_H_morning))
    # Count H desc array
    ds_TB.count_H_afternoon.attrs['long_name'] = 'Number of SMOS observations per pixel used in TB H polarization afternoon reconstruction'
    ds_TB.count_H_afternoon.attrs['standard_name'] = 'Number of H-pol afternoon obs per pixel'
    ds_TB.count_H_afternoon.attrs['units'] = 'Kelvins'
    ds_TB.count_H_afternoon.attrs['valid_min'] = '0'
    ds_TB.count_H_afternoon.attrs['valid_max'] = '9999'
    ds_TB.count_H_afternoon.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.count_H_afternoon),np.nanmax(ds_TB.count_H_afternoon))
    # Count V asc array
    ds_TB.count_V_morning.attrs['long_name'] = 'Number of SMOS observations per pixel used in TB V polarization morning reconstruction'
    ds_TB.count_V_morning.attrs['standard_name'] = 'Number of V-pol morning obs per pixel'
    ds_TB.count_V_morning.attrs['units'] = 'Kelvins'
    ds_TB.count_V_morning.attrs['valid_min'] = '0'
    ds_TB.count_V_morning.attrs['valid_max'] = '9999'
    ds_TB.count_V_morning.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.count_V_morning),np.nanmax(ds_TB.count_V_morning))
    # Count V desc array
    ds_TB.count_V_afternoon.attrs['long_name'] = 'Number of SMOS observations per pixel used in TB V polarization afternoon reconstruction'
    ds_TB.count_V_afternoon.attrs['standard_name'] = 'Number of V-pol afternoon obs per pixel'
    ds_TB.count_V_afternoon.attrs['units'] = 'Kelvins'
    ds_TB.count_V_afternoon.attrs['valid_min'] = '0'
    ds_TB.count_V_afternoon.attrs['valid_max'] = '9999'
    ds_TB.count_V_afternoon.attrs['actual_range'] = '[%.1f  %.1f]'%(np.nanmin(ds_TB.count_V_afternoon),np.nanmax(ds_TB.count_V_afternoon))
    # DGG multilayer fit parameters if multi-incidence
    if single_incidence==False:
        ds_TB.multilayer_fit_params.attrs['long_name'] = 'Optimal parameters and residual of the multilayer fit for incidence correction'
        ds_TB.multilayer_fit_params.attrs['standard_name'] = 'Optimal parameters and residual'
        ds_TB.multilayer_fit_params.attrs['structure'] = 'Columns 0 and 1: coordinates (X and Y) of the DGG. \nColumns 2 to 7: eps_1, eps_2, T, alpha, beta and residual of the fit for morning tracks. \nColumns 8 to 13: eps_1, eps_2, T, alpha, beta and residual of the fit for afternoon tracks.'
        ds_TB.DGG_ID.attrs['long_name'] = 'DGG ID from the L1C dataset'
        ds_TB.DGG_ID.attrs['standard_name'] = 'DGG_ID'
        
    # General attributes
    ds_TB.attrs['title'] = 'SMOS L-band enhanced resolution brightness temperatures based on rSIR image reconstruction technique'
    ds_TB.attrs['institution'] = "Institut des Géosciences de l'Environnement (CNRS/UGA) & Centre National d'Etudes Spatiales (CNES)"
    ds_TB.attrs['source'] = 'SMOS L1C Full polarization reconstructed TB swath over land (MIR_SCLF1C product)'
    ds_TB.attrs['product_version'] = 'v1.0'
    ds_TB.attrs['creation_date'] = '2024-09-02'
    ds_TB.attrs['created_by'] = 'Pierre Zeiger, IGE'
    ds_TB.attrs['contact'] = 'pierre.zeiger@univ-grenoble-alpes.fr'
    ds_TB.attrs['CRS'] = 'WGS 84 / NSIDC Sea Ice Polar Stereographic South'
    ds_TB.attrs['EPSG'] = '3976'
    ds_TB.attrs['spatial_sampling'] = '12500'
    ds_TB.attrs['spatial_units'] = 'meters'
    ds_TB.attrs['pixel_area'] = '156.25 km²'
    ds_TB.attrs['time_coverage_start'] = '2010-01-12'
    ds_TB.attrs['time_coverage_end'] = '2023-07-31'
    ds_TB.attrs['time_coverage_res'] = 'twice-daily'
    ds_TB.attrs['license'] = 'Restricted to the identified user'
    ds_TB.attrs['platform'] = 'SMOS'
    ds_TB.attrs['instrument'] = 'MIRAS'
    ds_TB.attrs['no_data_value'] = 'nan'
    ds_TB.attrs['num_obs_SMOS_init'] = '%s'%lens_data[0]
    ds_TB.attrs['num_obs_SMOS_flagged'] = '%s'%lens_data[1]
    ds_TB.attrs['num_obs_SMOS_flagged_morning'] = '%s'%len_asc_init
    ds_TB.attrs['num_obs_SMOS_flagged_afternoon'] = '%s'%len_desc_init
    ds_TB.attrs['num_DGG_SMOS'] = '%s'%len_DGG
    ds_TB.attrs['num_obs_rSIR_H_morning'] = '%s'%len_H_asc_fin
    ds_TB.attrs['num_obs_rSIR_H_afternoon'] = '%s'%len_H_desc_fin
    ds_TB.attrs['num_obs_rSIR_V_morning'] = '%s'%len_V_asc_fin
    ds_TB.attrs['num_obs_rSIR_V_afternoon'] = '%s'%len_V_desc_fin
    
    return ds_TB
    
    
    
def edit_tb_attributes(ds_TB, name):   # Change attributes dtype for matching the CF conventions (https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.pdf)
    # Remove all attributes
    attrs = [a for a in ds_TB.attrs]
    for attr in attrs:
        del ds_TB.attrs[attr] 
    
    # General attributes
    ds_TB.attrs['title'] = 'SMOS L-band enhanced resolution brightness temperatures based on rSIR image reconstruction technique'
    ds_TB.attrs['institution'] = "Institut des Géosciences de l'Environnement (CNRS/UGA) & Centre National d'Etudes Spatiales (CNES)"
    ds_TB.attrs['source'] = 'SMOS L1C Full polarization reconstructed TB swath over land (MIR_SCLF1C product), v724'
    ds_TB.attrs['product_version'] = 'v1.0'
    ds_TB.attrs['creation_date'] = '2024-09-02'
    ds_TB.attrs['created_by'] = 'Pierre Zeiger, IGE'
    ds_TB.attrs['contact'] = 'pierre.zeiger@univ-grenoble-alpes.fr'
    
    if name == 'Antartica':
        ds_TB.attrs['CRS'] = 'WGS 84 / NSIDC Sea Ice Polar Stereographic South'
        ds_TB.attrs['EPSG'] = 3976
    elif name == 'Greenland':
        ds_TB.attrs['CRS'] = 'WGS 84 / NSIDC Sea Ice Polar Stereographic North'
        ds_TB.attrs['EPSG'] = 3413
    ds_TB.attrs['spatial_sampling'] = 12500
    ds_TB.attrs['spatial_units'] = 'meter'
    ds_TB.attrs['pixel_area'] = 156.25
    ds_TB.attrs['pixel_area_unit'] = 'km²'
    ds_TB.attrs['time_coverage_start'] = '2010-04-12'
    ds_TB.attrs['time_coverage_end'] = '2024-03-31'
    ds_TB.attrs['time_coverage_res'] = 'twice-daily'
    ds_TB.attrs['license'] = 'Restricted to the identified user'
    ds_TB.attrs['platform'] = 'SMOS'
    ds_TB.attrs['instrument'] = 'MIRAS'
    ds_TB.attrs['no_data_value'] = np.nan
    
    # Data Array attributes
    for var in ds_TB.variables: 	# Change dtype to float
        if var in ['x', 'y', 'TB_H_morning', 'TB_H_afternoon', 'TB_V_morning', 'TB_V_afternoon']:
            ds_TB[var].attrs['long_name'] = ''.join(ds_TB[var].attrs['long_name'].split('_'))
            ds_TB[var].attrs['standard_name'] = ds_TB[var].attrs['standard_name']
            if var in ['x', 'y']:
                ds_TB[var].attrs['units'] = 'meter'
            else:
                ds_TB[var].attrs['units'] = 'Kelvin'
            ds_TB[var].attrs['valid_min'] = float(ds_TB[var].attrs['valid_min'])
            ds_TB[var].attrs['valid_max'] = float(ds_TB[var].attrs['valid_max'])
            ds_TB[var].attrs['actual_range'] = [float(ds_TB[var].attrs['actual_range'].split('  ')[0][1:]), float(ds_TB[var].attrs['actual_range'].split('  ')[1][:-1])]
            
        if var in ['iterations']:	# Change dtype to int
            ds_TB[var].attrs['long_name'] = ''.join(ds_TB[var].attrs['long_name'].split('_'))
            ds_TB[var].attrs['standard_name'] = ds_TB[var].attrs['standard_name']
            ds_TB[var].attrs['units'] = 1
            ds_TB[var].attrs['valid_min'] = int(ds_TB[var].attrs['valid_min'])
            ds_TB[var].attrs['valid_max'] = int(ds_TB[var].attrs['valid_max'])
            ds_TB[var].attrs['actual_range'] = [int(ds_TB[var].attrs['actual_range'].split('  ')[0][1:]), int(ds_TB[var].attrs['actual_range'].split('  ')[1][:-1])]
                        
        return ds_TB   #.drop_vars(['TB_H_acc_morning', 'TB_H_acc_afternoon', 'TB_V_acc_morning', 'TB_V_acc_afternoon', 'count_H_morning', 'count_V_morning', 'count_H_afternoon', 'count_V_afternoon'])
           
    
def edit_melt_attributes(ds_melt, name):   # Change attributes dtype for matching the CF conventions (https://cfconventions.org/Data/cf-conventions/cf-conventions-1.10/cf-conventions.pdf)
    # General attributes
    ds_melt.attrs['title'] = 'SMOS L-band enhanced resolution melt (wet/dry) status'
    ds_melt.attrs['institution'] = "Institut des Géosciences de l'Environnement (CNRS/UGA) & Centre National d'Etudes Spatiales (CNES)"
    ds_melt.attrs['source'] = 'SMOS rSIR-enhanced 40° incidence brightness temperatures, v1.0'
    ds_melt.attrs['product_version'] = 'v1.0'
    ds_melt.attrs['creation_date'] = '2024-09-02'
    ds_melt.attrs['created_by'] = 'Pierre Zeiger, IGE'
    ds_melt.attrs['contact'] = 'pierre.zeiger@univ-grenoble-alpes.fr'
    
    if name == 'Antartica':
        ds_melt.attrs['CRS'] = 'WGS 84 / NSIDC Sea Ice Polar Stereographic South'
        ds_melt.attrs['EPSG'] = 3976
    elif name == 'Greenland':
        ds_melt.attrs['CRS'] = 'WGS 84 / NSIDC Sea Ice Polar Stereographic North'
        ds_melt.attrs['EPSG'] = 3413
    ds_melt.attrs['spatial_sampling'] = 12500
    ds_melt.attrs['spatial_units'] = 'meter'
    ds_melt.attrs['pixel_area'] = 156.25
    ds_melt.attrs['pixel_area_unit'] = 'km²'
    ds_melt.attrs['time_coverage_start'] = '2010-04-12'
    ds_melt.attrs['time_coverage_end'] = '2024-03-31'
    ds_melt.attrs['time_coverage_res'] = 'twice-daily'
    ds_melt.attrs['license'] = 'Restricted to the identified user'
    ds_melt.attrs['platform'] = 'SMOS'
    ds_melt.attrs['instrument'] = 'MIRAS'
    ds_melt.attrs['no_data_value'] = np.nan
    
    # X array
    if name == 'Antartica':
        ds_melt.x.attrs['long_name'] = 'X coordinate NSIDC Polar Stereographic South'
    elif name == 'Greenland':
        ds_melt.x.attrs['long_name'] = 'X coordinate NSIDC Polar Stereographic North'
    ds_melt.x.attrs['standard_name'] = 'X'
    ds_melt.x.attrs['units'] = 'meter'
    ds_melt.x.attrs['valid_min'] = -3289335.29
    ds_melt.x.attrs['valid_max'] = 3289335.29
    ds_melt.x.attrs['actual_range'] = [np.nanmin(ds_melt.x.values), np.nanmax(ds_melt.x.values)]
    
    # Y array
    if name == 'Antartica':
        ds_melt.y.attrs['long_name'] = 'Y coordinate NSIDC Polar Stereographic South'
    elif name == 'Greenland':
        ds_melt.y.attrs['long_name'] = 'Y coordinate NSIDC Polar Stereographic North'
    ds_melt.y.attrs['standard_name'] = 'Y'
    ds_melt.y.attrs['units'] = 'meter'
    ds_melt.y.attrs['valid_min'] = -3323160.27
    ds_melt.y.attrs['valid_max'] = 3323160.27
    ds_melt.y.attrs['actual_range'] = [np.nanmin(ds_melt.y.values), np.nanmax(ds_melt.y.values)]
    
    # snow_status_wet_dry_smos_morning array
    ds_melt.snow_status_wet_dry_smos_morning.attrs['long_name'] = 'Snow status SMOS morning tracks wet or dry'
    ds_melt.snow_status_wet_dry_smos_morning.attrs['standard_name'] = 'snow_status_wet_dry_smos_morning'
    ds_melt.snow_status_wet_dry_smos_morning.attrs['units'] = 1
    ds_melt.snow_status_wet_dry_smos_morning.attrs['valid_min'] = 0
    ds_melt.snow_status_wet_dry_smos_morning.attrs['valid_max'] = 1
    ds_melt.snow_status_wet_dry_smos_morning.attrs['actual_range'] = [np.nanmin(ds_melt.snow_status_wet_dry_smos_morning.values), np.nanmax(ds_melt.snow_status_wet_dry_smos_morning.values)]
    
    # snow_status_wet_dry_smos_afternoon array
    ds_melt.snow_status_wet_dry_smos_afternoon.attrs['long_name'] = 'Snow status SMOS afternoon tracks wet or dry'
    ds_melt.snow_status_wet_dry_smos_afternoon.attrs['standard_name'] = 'snow_status_wet_dry_smos_afternoon'
    ds_melt.snow_status_wet_dry_smos_afternoon.attrs['units'] = 1
    ds_melt.snow_status_wet_dry_smos_afternoon.attrs['valid_min'] = 0
    ds_melt.snow_status_wet_dry_smos_afternoon.attrs['valid_max'] = 1
    ds_melt.snow_status_wet_dry_smos_afternoon.attrs['actual_range'] = [np.nanmin(ds_melt.snow_status_wet_dry_smos_afternoon.values), np.nanmax(ds_melt.snow_status_wet_dry_smos_afternoon.values)]
    
    return ds_melt
    
























