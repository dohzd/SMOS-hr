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


if __name__ == "__main__":
    # No warnings
    warning_off = True
    if warning_off:
        warnings.filterwarnings("ignore", category=RuntimeWarning)
    np.set_printoptions(suppress=True)  # Put off scientific notation

    now = time.time()

    root_dir = os.environ["ROOT_DIR"]
    dl_dir = os.environ["DL_DIR"]
    area = os.environ["AREA"]
    grid_choice = os.environ["GRID"]
    assert (area=='Greenland') | (area=='Antarctica'), "$NAME variable should be one of: 'Greenland', 'Antarctica'"
    grid = xr.open_dataset(f"{root_dir}/mask/{area}_grid_{grid_choice}km_lonlat.nc")
    # Reconstruct
    print("Starting reconstruction")
    ds_TB = extract_mat.reconstruct_from_mat(dl_dir, grid, subset=None)

    print('Finished. Total time ellapsed: %.2f hours'%((time.time()-now)/3600))





