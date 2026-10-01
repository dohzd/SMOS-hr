# -*- coding: utf-8 -*-
"""
Created on Thu May 11 15:31:33 2023

@author: Pierre
"""

import numpy as np
import xarray as xr
import time 
from tqdm import tqdm
from scipy.spatial import KDTree


def compute_MRF_list(x,y,x0,y0,rotation,semi_minor,semi_major,pixel_list):
    # Compute MRF from distance to the measurement center and azimuth angle
    # Returns two lists (the first is the list of pixels in the ellispe, the second is the MRF i.e. the weighted values)
    # Returning lists instead of 2D arrays saves memory up to a factor 10000 	
    term1 = (((x - x0) * np.cos(rotation*np.pi/180-np.pi/2) - (y - y0) * np.sin(rotation*np.pi/180-np.pi/2))**2 / semi_major**2)
    term2 = (((x - x0) * np.sin(rotation*np.pi/180-np.pi/2) + (y - y0) * np.cos(rotation*np.pi/180-np.pi/2))**2 / semi_minor**2)
    ellipse = (term1 + term2) <= 1
    MRF = np.asarray(np.log(1/2) * np.exp(4*term1[ellipse] + 4*term2[ellipse]))
    return pixel_list[ellipse].astype(int),MRF#[ellipse]
    

def extract_MRFs(x_ref,y_ref,X,Y,a,b,azimuth,landseamask=None,resolution=12500):   # Compute MRFs
    pixel_list = np.arange(len(x_ref)*len(y_ref)).astype(float)	# pixel's number
    xx = np.repeat(x_ref[np.newaxis,:],len(y_ref),axis=0).flatten() + resolution/2    # 6.25 km shift to move (x,y) from left corner to pixels centers
    yy = np.repeat(y_ref[:,np.newaxis],len(x_ref),axis=1).flatten() - resolution/2    # 6.25 km shift to move (x,y) from left corner to pixels centers
    MRFs_list = np.asarray([compute_MRF_list(xx,yy,X[i],Y[i],azimuth[i],b[i],a[i],pixel_list) for i in range(len(X))],dtype=object)
    return MRFs_list 


def AVE_iter(BT_H,BT_V,BT_H_acc,BT_V_acc,mask_H_acc,mask_V_acc,gain):
    denom_H = np.asarray([10**(gain/20)*BT_H,np.repeat(0,len(gain))])[np.asarray([mask_H_acc, not mask_H_acc])]
    numer_H = np.asarray([10**(gain/20),np.repeat(0,len(gain))])[np.asarray([mask_H_acc, not mask_H_acc])]
    denom_H_acc = np.asarray([10**(gain/20)*BT_H_acc,np.repeat(0,len(gain))])[np.asarray([mask_H_acc, not mask_H_acc])]
    denom_V = np.asarray([10**(gain/20)*BT_V,np.repeat(0,len(gain))])[np.asarray([mask_V_acc, not mask_V_acc])]
    numer_V = np.asarray([10**(gain/20),np.repeat(0,len(gain))])[np.asarray([mask_V_acc, not mask_V_acc])]
    denom_V_acc = np.asarray([10**(gain/20)*BT_V_acc,np.repeat(0,len(gain))])[np.asarray([mask_V_acc, not mask_V_acc])]
    count_H = np.asarray([np.ones(len(gain)),np.zeros(len(gain))])[np.asarray([mask_H_acc, not mask_H_acc])]
    count_V = np.asarray([np.ones(len(gain)),np.zeros(len(gain))])[np.asarray([mask_V_acc, not mask_V_acc])]
    return np.asarray([denom_H,numer_H,denom_H_acc,denom_V,numer_V,denom_V_acc,count_H,count_V]).squeeze()
    

def AVE_direct(x_ref,y_ref,X,Y,MRFs,BT_H,BT_V,BT_H_acc,BT_V_acc,mask_H_acc=None,mask_V_acc=None,dist_pix=15):     # 
    if len(MRFs)>0:    
        if mask_H_acc is None:
            mask_H_acc = np.ones(len(MRFs)).astype(bool)
        if mask_V_acc is None:
            mask_V_acc = np.ones(len(MRFs)).astype(bool) 
        mask_pix = mask_pixel_dist(x_ref,y_ref,X[mask_V_acc],Y[mask_V_acc],thres=dist_pix,resolution=12500) 	# Mask pixels based on V-pol positions (i.e. the most conservative as V_acc > H_acc)
        data_vect = np.repeat(np.zeros(len(x_ref)*len(y_ref))[np.newaxis,:],8,axis=0)
        for i,m in enumerate(MRFs):   # Compute AVE for H pol with mask based on accuracy
            data_vect[:,m[0].astype(int)] = data_vect[:,m[0].astype(int)] + AVE_iter(BT_H[i],BT_V[i],BT_H_acc[i],BT_V_acc[i],mask_H_acc[i],mask_V_acc[i],m[1])
        data_vect[np.repeat((~mask_pix)[np.newaxis,:],8,axis=0)] = np.nan
        AVE_H,AVE_H_acc,AVE_V,AVE_V_acc,count_H,count_V = [data_vect[0]/data_vect[1],data_vect[2]/data_vect[1],
            data_vect[3]/data_vect[4],data_vect[5]/data_vect[4],data_vect[6],data_vect[7]]
            
    else:
        AVE_H,AVE_H_acc,AVE_V,AVE_V_acc,count_H,count_V = np.zeros((6,len(y_ref)*len(x_ref)))*np.nan
    return AVE_H,AVE_V,AVE_H_acc,AVE_V_acc,count_H,count_V#,mask_pix_H,mask_pix_V



def mask_pixel_dist(x_ref,y_ref,X,Y,thres=10,resolution=12500):
    time0 = time.time()
    xx = (np.repeat(x_ref[np.newaxis,:],len(y_ref),axis=0) + resolution/2).flatten()   # 6.25 km shift to move (x,y) from left corner to pixels centers
    yy = (np.repeat(y_ref[:,np.newaxis],len(x_ref),axis=1) - resolution/2).flatten()   # 6.25 km shift to move (x,y) from left corner to pixels centers
    tree = KDTree(np.transpose([X,Y]),leafsize=100)
    dist = tree.query(np.transpose([xx,yy]))[0]
    print("Pixel mask with a %s km threshold: %.2f hours ellapsed"%(thres,(time.time()-time0)/3600))
    return np.asarray(dist)<=thres*1000
    
    
def update_term(Tb_guess,forward,scale,gain):
    if scale>=1:
        update = 1/(1/(2*forward) * (1-1/scale) + 1/(Tb_guess*scale))
    else:
        update = (1/2) * forward * (1-scale) + Tb_guess*scale
    return np.asarray([10**(gain/20),10**(gain/20)*update]).squeeze()
    

def rSIR_iter(MRFs,BT,AVE,max_iter=5,mask_acc=None):
    Tb_guess = AVE
    if mask_acc is not None:
        BT = BT[mask_acc]
        MRFs = MRFs[mask_acc]
    n_iter = 1
    Tb_estimates = [Tb_guess]
    t1 = time.time()
    while (n_iter<=max_iter):# & (conv_std>0.001):
        mask_nan = np.isnan(Tb_guess)   
        forward = np.asarray([(1/np.nansum(10**(MRFs[i][1]/20)[np.asarray([~mask_nan[m] for m in MRFs[i][0]])])) * 
                              np.nansum(10**(MRFs[i][1]/20) * [Tb_guess[m] for m in MRFs[i][0]]) for i in range(len(MRFs))])
        forward[forward==0] = np.nan
        scale = np.sqrt(BT/forward)
        sum_MRF = np.repeat(np.zeros(len(Tb_guess))[np.newaxis,:],2,axis=0)            
        for i,m in enumerate(MRFs):      
            sum_MRF[:,m[0]] = sum_MRF[:,m[0]] + update_term(Tb_guess[m[0]],forward[i],scale[i],m[1])
            
        Tb_guess_update = sum_MRF[1]/sum_MRF[0]
        Tb_guess_update[np.isnan(Tb_guess)] = np.nan
        Tb_estimates.append(Tb_guess_update)
        Tb_guess = Tb_guess_update
        n_iter = n_iter+1
    dt = time.time()-t1
    print('rSIR finished: %.0f seconds per iteration'%(dt/15))
    Tb_estimates = np.asarray(Tb_estimates)   
    
    return Tb_estimates
    




    

