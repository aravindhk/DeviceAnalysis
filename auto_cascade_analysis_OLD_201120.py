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
from scipy.signal import find_peaks  # For finding local maxima in gm/dgm

# Changing plot properties
plt.style.use(['default', 'seaborn-bright'])
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['xtick.labelsize'] = 9
mpl.rcParams['ytick.labelsize'] = 9
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 10
mpl.rcParams['lines.markersize'] = 3
mpl.rcParams['figure.subplot.wspace'] = 0.3  # Horizontal spacing in subplots
# mpl.rcParams['figure.subplot.hspace'] = 0.3  # Vertical spacing in subplots
mpl.rcParams['figure.subplot.hspace'] = 0.4  # Vertical spacing in subplots
mpl.rcParams['figure.figsize'] = 12, 6  # Figure size
mpl.rcParams['figure.constrained_layout.h_pad'] = 0.04167  # Padding around axes objects
mpl.rcParams['figure.constrained_layout.w_pad'] = 0.04167
mpl.rcParams['legend.fontsize'] = 'small'
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['legend.borderpad'] = 0


import pandas as pd
import xlrd
import csv
import numpy as np
import seaborn as sns
import statsmodels.api as sm
from pathlib import Path
from device_analysis_classes import *

# Close all plots from previous run(s)
plt.close('all')

# Choose plotting colors here
pal = sns.color_palette("bright")
custom_colors = pal.as_hex()
custom_colors = custom_colors[:5]+custom_colors[6:7]+custom_colors[-1:]
# custom_colors = custom_colors[3:]
sns.set_palette(custom_colors)
# sns.palplot(sns.color_palette("bright", 8))
# custom_colors = custom_colors[:1]+custom_colors[2:5]+custom_colors[6:7]+custom_colors[-1:]
# custom_colors = ["DarkBlue", "Crimson", "DarkGreen", "DarkOrchid", "Orange", "DeepPink"]
# You can use any HTML colors in the above line. For HTML colors, see this: https://www.w3schools.com/colors/colors_names.asp


# Setting current directory and target directory
user_folder = "F:\\Google Drive\\Research\\Projects"
# user_folder = "C:\\Users\\aravi\\Google Drive\\Research\\Projects"

##### Selecting the folder ####################################################
# dir_path = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-03-04-MF1"; Vg_idvd_plt = 40; Vds_low = 5;
# dir_path = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-07-11-MF1"; Vg_idvd_plt = 40; Vds_low = 5;
# dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-01-IB13"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 0;
# dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-01-IB14"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 0;
# dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-06-IB13-source-ground"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-01-IB13-N2-anneal"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-03-IB13-N2-anneal"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# dir_path = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-02-MH3"; Vg_idvd_plt = 30; Vds_low = 0.5; isbipolar_idvg = 1;
# dir_path = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-13-MJ2"; Vg_idvd_plt = 30; Vds_low = 0.5; isbipolar_idvg = 1;
# dir_path = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-13-MJ3"; Vg_idvd_plt = 30; Vds_low = 0.5; isbipolar_idvg = 1;
dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-19-IB21A"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
###############################################################################

target_dir = user_folder + dir_path
os.chdir(target_dir) #Move to target directory

#Reading file parameters from input files
#Reading channel parameters - how many channels, what channel names, etc.
channel_params = Path("channel_parameters.xlsx") #Channel parameters file
channel_params = xlrd.open_workbook(channel_params,logfile=open(os.devnull,'w'))
channel_params = np.asarray(pd.read_excel(channel_params,engine='xlrd')) #Reading Id-Vg data
channel_length = channel_params[:,0] #Channel lengths (string)
channel_L = channel_params[:,1] #Channel lengths (um)
channel_index = np.arange(channel_L.shape[0]) #Channel lengths index

# Constants
Cox_30nm = 116 #30nm SiO2 capacitance (Kirby ACS Nano) (nF/cm2)

#Reading device parameters
device_params = Path("device_parameters_auto.xlsx") #Device parameters file - what sweeps, column indices, etc.
device_params = xlrd.open_workbook(device_params,logfile=open(os.devnull,'w'))
device_params = np.asarray(pd.read_excel(device_params,engine='xlrd')) #Reading Id-Vg data
chip = device_params[:,2] #List of chips (needed?)
device = device_params[:,2] #List of devices
file_prefix = device_params[:,3] #List of file prefixes
isfolder_organized = device_params[:,4] #Whether the device data is organized into folders

W = device_params[0,14] #Channel width (um)
tox = device_params[0,15] #Thickness of the gate oxide (nm)
Cox = Cox_30nm*30/tox #100nm SiO2 capacitance (nF/cm2)

device_index = np.arange(device.shape[0]) #Device index


# Main section of code
os.chdir(target_dir)
tlm_set = [None]*device_index.shape[0] #Creating None array to hold TLM objects

device_count = 0

fig, axs = plt.subplots(2,3)  
fig2, axs2 = plt.subplots(2,3)  
fig3, axs3 = plt.subplots(2,3)  
fig4, axs4 = plt.subplots(2,3)  
axs_reshape = axs.reshape(-1)
axs2_reshape = axs2.reshape(-1)
axs3_reshape = axs3.reshape(-1)
axs4_reshape = axs4.reshape(-1)
# fig, ax = plt.subplots(2,2)
# ax = ax.reshape(-1)

for i in device_index:    
    # print(device_count)
    if pd.isnull(file_prefix[i]): #Empty file prefix is read as NaN from excel
        file_prefix[i] = ""
    if isfolder_organized[i]: #If organized into folders, cd into folder
        os.chdir(device[i])        
    col_index = device_params[np.ix_([i],[5,6,7,8,9,10,19])] #Column indices for various values
    col_index = col_index[0] #Converting 2-D array to 1-D
    #1,2,3 - Ids, Vgs, Vds for Id-Vg
    #4,5,6 - Ids, Vgs, Vds for Id-Vd    
    num_vd_idvg = device_params[i, 16]  # Number of Vd steps in Id-Vg sweep
    num_var_idvg = device_params[i, 17]  # Number of column variables in a single Id-Vg sweep
    # num_vg_idvd = device_params[i, 11]  # Number of Vg steps in Id-Vd sweep
    # isbipolar_idvd = device_params[i, 12]  # Whether Id-Vd sweep is bipolar
    # num_var_idvd = device_params[i, 13]  # Number of column variables in a single Id-Vd sweep

    W = device_params[i, 14]  # Channel width (um)
    tox = device_params[i, 15]  # Thickness of the gate oxide (nm)
    Cox = Cox_30nm * 30 / tox  # 100nm SiO2 capacitance (nF/cm2)
    input_data = [None] * channel_length.shape[0]  # Creating dummy array for input

    # channel_select = list("C111111")
    channel_select = list(device_params[i, 18]) #Encoding to include/exclude channels

    channel_count = 0    
    for j in channel_index:                    
        filename_g = file_prefix[i] + channel_length[j] + "-G-Vds" + str(Vds_low) + ".csv" #Filename of id-vg files
#        filename_d = file_prefix[i] + channel_length[j] + "-D.xls" #Filename of id-vd file
        #Reading Id-Vg files
        my_file_g = Path(filename_g) #Creating file object        
#        my_file_d = Path(filename_d) #Creating file object

        if my_file_g.exists() and int(channel_select[j+1]): #Checking if IdVg and IdVd files exist
#            wb_idvg = xlrd.open_workbook(filename_g,logfile=open(os.devnull,'w'))                     
            # print(my_file_g)
            csv_file = open(my_file_g)            
            idvg_data = np.asarray(pd.read_csv(my_file_g)) #Reading Id-Vg data            
            input_data[channel_count] = InputData(idvg_data, channel_length, channel_L, j) #Assigning the read data to InputData object
            channel_count += 1                 

    if channel_count:        
        params = Parameters(col_index, num_var_idvg, num_vd_idvg, W, Cox, isbipolar_idvg = isbipolar_idvg, Vds_param = Vds_low)
        tlm_set[device_count] = TLM(input_data[:channel_count],channel_count, params) #Creating TLM object
        tlm_set[device_count].tlm_calc(params) #Performing TLM calculations            
        if tlm_set[device_count].R_squared>0.969:
            print(tlm_set[device_count].R_squared)
            # print(filename_g)
            # mpl.rcParams['figure.figsize'] = 12, 6
            # fig, ((ax1, ax2, ax5),(ax3, ax4, ax6)) = plt.subplots(2,3)  
            # fig, ax1 = plt.subplots(1,1)  
            
            # fig_2, ax_2 = plt.subplots(1,1)
            # fig_2.set_size_inches(4, 3, forward=True)
            # tlm_set[i].plot_IdVg(ax_2, fig_2, lg = 'log', vd_step = 'low')  # Plotting Id-Vg for all channel lengths in log scale         
            # fig_2.suptitle(device[i], fontsize=14)
            # ax_2.grid()            
            tlm_set[device_count].plot_IdVg(axs_reshape[device_count],fig,'log') #Plotting Id-Vg for all channel lengths in log scale
            # tlm_set[device_count].plot_IdVg(ax[0],fig,'log') #Plotting Id-Vg for all channel lengths in log scale    
            axs_reshape[device_count].set_title(device[i])
            tlm_set[device_count].plot_Rc(axs2_reshape[device_count],fig2,'n2D') #Plotting Rc vs n2Ds
            # tlm_set[device_count].plot_Rc(ax[1],fig,'n2D') #Plotting Rc vs n2Ds
            axs2_reshape[device_count].set_title(device[i])
            tlm_set[device_count].plot_mu_TLM(axs3_reshape[device_count],fig3,'n2D') #Plotting mu_TLM vs n2D
            # tlm_set[device_count].plot_mu_TLM(ax[2],fig,'n2D') #Plotting mu_TLM vs n2D
            axs3_reshape[device_count].set_title(device[i])
            tlm_set[device_count].plot_TLM_fit(axs4_reshape[device_count],fig4) #Plotting TLM fits
            # tlm_set[device_count].plot_TLM_fit(ax[3],fig) #Plotting TLM fits
            axs4_reshape[device_count].set_title(device[i])
            # tlm_set[device_count].plot_Rsh(ax5, fig, 'n2D')  # Plotting Rsh vs n2D        
            # tlm_set[device_count].plot_IdVg_Vov(ax6,fig,'log') #Plotting Id-Vg for all channel lengths in log scale    
            # print(np.amax(tlm_set[device_count].Id_adjusted[0,:])*1e6)
            device_count = device_count + 1            
            # fig.savefig("summary-" + str(tlm_set[i].R_squared) + ".svg", \
            #                 transparent=True, bbox_inches='tight', pad_inches=0.1)
            # fig_2.savefig("summary_idvg-" + str(i) + ".svg", transparent=True,\
            #                             bbox_inches='tight', pad_inches=0.1)    
    
    if isfolder_organized[i]:
        os.chdir(target_dir)

fig.savefig("summary-idvg.svg", transparent=True, bbox_inches='tight', pad_inches=0.1)
fig2.savefig("summary-rc.svg", transparent=True, bbox_inches='tight', pad_inches=0.1)
fig3.savefig("summary-mutlm.svg", transparent=True, bbox_inches='tight', pad_inches=0.1)
fig4.savefig("summary-tlmfit.svg", transparent=True, bbox_inches='tight', pad_inches=0.1)
plt.show()            