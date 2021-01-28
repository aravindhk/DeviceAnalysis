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

# Setting current directory and target directory
user_folder = "F:\\Google Drive\\Research\\Projects"
# user_folder = "C:\\Users\\aravi\\Google Drive\\Research\\Projects"
##### Selecting the folder ####################################################
# target_dir = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-03-04-MF1\\"
# target_dir = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-07-11-MF1\\"; W = 10; tox = 100; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-07-22-IB3\\Old\\"
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-01-IB13\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-01-IB14\\"; W = 2; tox = 80;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-06-IB13-source-ground\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-01-IB13-N2-anneal\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-03-IB13-N2-anneal\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-20-IB18-N2-purge\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-19-IB21A\\"; W = 2; tox = 50;
target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-24-IB21A\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-19-IB21B\\"; W = 2; tox = 50;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-02-MH3\\"; W = 1; tox = 95;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-13-MJ2\\"; W = 1; tox = 95;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-13-MJ3\\"; W = 1; tox = 95;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-23-MJ1\\"; W = 1; tox = 95;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-29-MJ5\\"; W = 1; tox = 95;
###############################################################################

# num_rows = 13
num_rows = 15
# num_cols = 10
num_cols = 9

num_device = 0
device_number = []
device = []

path = user_folder + target_dir

file_prefix = []
folder_organized = 0
Ids_IdVg = 2
Vgs_IdVg = 3
Vds_IdVg = 6
Igs_IdVg = 0
Ids_IdVd = 1
Vgs_IdVd = 6
Vds_IdVd = 2
num_Vg_IdVd = 4
bipolar_IdVd = 1
num_var_IdVd = 7
# num_Vd_IdVg = 1
# num_var_IdVg = 4
channel_select = "C111111"

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
        "channel_width(um)":W,"oxide_thickness(nm)":tox,"num_Vd_IdVg":num_Vd_IdVg,
        "num_var_IdVg":num_var_IdVg,"Channel_select":channel_select,"Igs_IdVg":Igs_IdVg}
df = pd.DataFrame(data)
df.to_excel(path+"device_parameters_auto.xlsx",index=False)

channel_labels = ["CH1", "CH2", "CH3", "CH4", "CH5", "CH6"]
channel_lengths = [0.1, 0.2, 0.3, 0.5, 0.7, 1]
channel_data = {"Channel Lengths":channel_labels, "Channel L":channel_lengths}
df = pd.DataFrame(channel_data)
df.to_excel(path+"channel_parameters.xlsx",index=False)