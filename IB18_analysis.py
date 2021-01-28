"""
Created on Thursday September 17 20:15:00 2020

@author: Aravindh Kumar
"""

"""
Code to analyse auto-cascade data
Currently only Id-Vg analysis is supported
"""

import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
import seaborn as sns
from pathlib import Path
from device_analysis_classes import *

# Changing plot properties
plt.style.use(['default', 'seaborn-bright'])
mpl = plot_properties(mpl)

# Close all plots from previous run(s)
plt.close('all')

# Choose plotting colors here
pal = sns.color_palette("bright")
custom_colors = pal.as_hex()
custom_colors = custom_colors[:5]+custom_colors[6:7]+custom_colors[-1:]
sns.set_palette(custom_colors)
# custom_colors = ["DarkBlue", "Crimson", "DarkGreen", "DarkOrchid", "Orange", "DeepPink"]
# You can use any HTML colors in the above line. For HTML colors, see this: https://www.w3schools.com/colors/colors_names.asp

# Setting current directory and target directory
user_folder = "F:\\Google Drive\\Research\\Projects"
# user_folder = "C:\\Users\\aravi\\Google Drive\\Research\\Projects"
##### Selecting the folder ####################################################
# dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-06-IB13-source-ground"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# dir_path = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-02-MH3"; Vg_idvd_plt = 30; Vds_low = 0.5; isbipolar_idvg = 1;
# dir_path = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-13-MJ2"; Vg_idvd_plt = 30; Vds_low = 0.5; isbipolar_idvg = 1;
# dir_path = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-13-MJ3"; Vg_idvd_plt = 30; Vds_low = 0.5; isbipolar_idvg = 1;
# dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-06-IB13-source-ground\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-20-IB18-N2-purge\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-24-IB21A\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# dir_path = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-23-MJ1"; Vg_idvd_plt = 30; Vds_low = 0.5; isbipolar_idvg = 1;
# dir_path = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-29-MJ5"; Vg_idvd_plt = 30; Vds_low = 1; isbipolar_idvg = 1;
# dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-06-IB13-source-ground"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# dir_path = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-07-11-MF1"; Vg_idvd_plt = 40; Vds_low = 5; isbipolar_idvg = 0;
###############################################################################
target_dir = user_folder + dir_path
os.chdir(target_dir) #Move to target directory

#Reading file parameters from input files
channel_params_file = "channel_parameters.xlsx"
device_params_file = "device_parameters_auto.xlsx"
channel_params, device_params = read_params(channel_params_file, device_params_file)

# Constants
Cox_30nm = 116 #30nm SiO2 capacitance (Kirby ACS Nano) (nF/cm2)

#Channel parameters
channel_length = channel_params[:,0] #Channel lengths (string)
channel_L = channel_params[:,1] #Channel lengths (um)

#Device parameters
device = device_params[:,2] #List of devices
file_prefix = device_params[:,3] #List of file prefixes
isfolder_organized = device_params[:,4] #Is data organized into folders?
W = device_params[0,14] #Channel width (um)
tox = device_params[0,15] #Thickness of the gate oxide (nm)
Cox = Cox_30nm*30/tox #100nm SiO2 capacitance (nF/cm2)

# Main section of code
os.chdir(target_dir)
tlm_set = [None]*device.shape[0] #Creating None array to hold TLM objects

device_count = 0

fig, ax = plt.subplots(2,3)
ax = ax.reshape(-1)
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()

gate_leak = np.zeros((15,9))

for i in np.arange(device.shape[0]):        
    if pd.isnull(file_prefix[i]): #Empty file prefix is read as NaN from excel
        file_prefix[i] = ""
    if isfolder_organized[i]: #If organized into folders, cd into folder
        os.chdir(device[i])        
    col_index = device_params[np.ix_([i],[5,6,7,8,9,10,19])] #Column indices for various values
    col_index = col_index[0] #Converting 2-D array to 1-D    
    num_vd_idvg = device_params[i, 16]  # Number of Vd steps in Id-Vg sweep
    num_var_idvg = device_params[i, 17]  # Number of column variables in a single Id-Vg sweep    
    W = device_params[i, 14]  # Channel width (um)
    tox = device_params[i, 15]  # Thickness of the gate oxide (nm)
    Cox = Cox_30nm * 30 / tox  # 100nm SiO2 capacitance (nF/cm2)
    input_data = [None] * channel_length.shape[0]  # Creating dummy array for input    
    channel_select = list(device_params[i, 18]) #Encoding to include/exclude channels

    channel_count = 0    
    for j in np.arange(channel_L.shape[0]):                            
        my_file_g = Path(file_prefix[i] + channel_length[j] + "-G-Vds" + str(Vds_low) + ".csv") #Filename of id-vg files
        # my_file_g = Path(file_prefix[i] + channel_length[j] + "-G.csv") #Filename of id-vg files
#        my_file_d = Path(file_prefix[i] + channel_length[j] + "-D.xls") #Filename of id-vd file        
        # print(my_file_g)
        if my_file_g.exists() and int(channel_select[j+1]): #Checking if IdVg and IdVd files exist
            # print(my_file_g)
            csv_file = open(my_file_g)            
            idvg_data = np.asarray(pd.read_csv(my_file_g)) #Reading Id-Vg data            
            if np.amax(idvg_data) > 1e101:
                pass
            else:
                input_data[channel_count] = InputData(idvg_data, channel_length, channel_L, j) #Assigning the read data to InputData object                        
                channel_count += 1

    if channel_count:            
        params = Parameters(col_index, num_var_idvg, num_vd_idvg, W, Cox, isbipolar_idvg = isbipolar_idvg, Vds_param = Vds_low)
        # print("device="+str(device_count))
        # print("channel="+str(channel_count))
        # print(my_file_g)
        tlm_set[device_count] = TLM(input_data[:channel_count],channel_count, params) #Creating TLM object        
        
        if tlm_set[device_count].gateLeak:
            device_string = device[i]
            r,c = device_string.replace('-', ' ').split(' ')
            gate_leak[int(r)-1, int(c)-1] = 1

        if 1:#tlm_set[device_count].R_squared>0.8:                                
            if tlm_set[device_count].gateLeak == 0:                
                tlm_set[device_count].tlm_calc(params) #Performing TLM calculations    
                if tlm_set[device_count].goodTLM and tlm_set[device_count].R_squared>0.8:
                    if tlm_set[device_count].count > 1:                        
                        mpl.rcParams['figure.figsize'] = 15, 3
                        fig_x, (ax1_x, ax2_x, ax3_x, ax4_x) = plt.subplots(1, 4) #More than 1 channel works
                        tlm_set[device_count].plot_IdVg(ax1_x, fig_x, lg = 'log', vd_step = 'hi')  # Plotting Id-Vg for all channel lengths in log scale
                        # tlm_set[device_count].plot_IdVd(ax2_x, fig_x, params = params, Vgs = Vg_idvd_plt)  # Plotting Id-Vd for all channel lengths @ Vgs=50V
                        # tlm_set[device_count].idvd[0].plot_IdVd(ax3_x, fig_x, params)  # Plotting Id-Vd for 1st channel length
                        tlm_set[device_count].plot_TLM_fit(ax2_x, fig_x)  # Plotting TLM fit for the last Vov
                        tlm_set[device_count].plot_Rc(ax3_x, fig_x, flag = 'n2D')  # Plotting Rc vs n2D
                        tlm_set[device_count].plot_mu_TLM(ax4_x, fig_x, flag = 'n2D')  # Plotting mu_eff vs n2D
                    else:                              
                        mpl.rcParams['figure.figsize'] = 5, 3
                        fig_x, ax1_x = plt.subplots(1, 1) #Only 1 channel works
                        tlm_set[device_count].plot_IdVg(ax1_x, fig_x, lg = 'log', vd_step = 'hi')  # Plotting Id-Vg for all channel lengths in log scale
                        # tlm_set[device_count].plot_IdVd(ax2_x, fig_x, params = params, Vgs = Vg_idvd_plt)  # Plotting Id-Vd for all channel lengths @ Vgs=50V
                        # tlm_set[device_count].idvd[0].plot_IdVd(ax3_x, fig_x, params)  # Plotting Id-Vd for 1st channel length                            

            for k in np.arange(channel_count):
                ax_index = np.where(channel_L == tlm_set[device_count].L[k])[0]        
                ax_index = ax_index[0]       

                if tlm_set[device_count].idvg[k,0].gateLeak == 0:                    
                    tlm_set[device_count].plot_IdVg_multVds(ax[ax_index], fig, params, 'log', channel = k, anno = 'yes')                      
                tlm_set[device_count].idvg[k,0].plot_Ig(ax2, fig2)  # Plotting Id-Vg for 100nm in log scale for smallest and largest vds           
            device_count = device_count + 1                         

    if isfolder_organized[i]:
        os.chdir(target_dir)
            

heat_map(ax3, gate_leak,colormap = 'seismic') #Heatmap of gate leakage
fig3.set_size_inches(6,6)
# fig3.savefig('gateleak_map.svg',transparent=True, bbox_inches='tight', pad_inches=0.1)

ax2.grid()
fig2.set_size_inches(8,6)
fig.savefig("summary-idvg.svg", transparent=True, bbox_inches='tight', pad_inches=0.1)
fig2.savefig("gate-leakage.svg", transparent=True, bbox_inches='tight', pad_inches=0.1)

plt.show()            