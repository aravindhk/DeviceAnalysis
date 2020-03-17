# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 04:16:18 2020

@author: Aravindh Kumar
"""

"""
Code to create device parameter file for auto-cascade data
"""

import numpy as np
import pandas as pd

num_rows = 8
num_cols = 8

num_device = 0
device_number = []
device = []
path = "\\Google Drive\\Research\\Projects\\InGaAs contacts\\Semi-Auto Cascade\\2020-03-04-MF1"
file_prefix = []
folder_organized = 0
Ids_IdVg = 2
Vgs_IdVg = 3
Vds_IdVg = 6
Ids_IdVd = 1
Vgs_IdVd = 6
Vds_IdVd = 2
num_Vg_IdVd = 4
bipolar_IdVd = 1
num_var_IdVd = 7
channel_width = 10
oxide_thickness = 100

for i in np.arange(num_rows):
    for j in np.arange(num_cols):
        num_device = num_device + 1
        device_number.append(num_device)
        device.append(str(i+1)+'-'+str(j+1))
        file_prefix.append(str(i+1)+'-'+str(j+1)+'-')
data = {"Device number":device_number,"Path":path,"Device":device,"File Prefix":file_prefix,
        "Folder Organized":folder_organized,"Ids_IdVg":Ids_IdVg, "Vgs_IdVg":Vgs_IdVg,
        "Vds_IdVg":Vds_IdVg,"Ids_IdVd":Ids_IdVd,"Vgs_IdVd":Vgs_IdVd,"Vds_IdVd":Vds_IdVd,
        "num_Vg_IdVd":num_Vg_IdVd,"bipolar_IdVd":bipolar_IdVd,"num_var_IdVd":num_var_IdVd,
        "channel_width(um)":channel_width,"oxide_thickness(nm)":oxide_thickness}
df = pd.DataFrame(data)
df.to_excel('device_parameters_auto.xlsx',index=False)