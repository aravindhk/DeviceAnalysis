# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 03:30:41 2020

@author: Aravindh Kumar
"""

"""
Code to analyse Cascade/Janis data
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
mpl.rcParams['figure.subplot.hspace'] = 0.3  # Vertical spacing in subplots
mpl.rcParams['figure.figsize'] = 12, 6  # Figure size
mpl.rcParams['figure.constrained_layout.h_pad'] = 0.04167  # Padding around axes objects
mpl.rcParams['figure.constrained_layout.w_pad'] = 0.04167
mpl.rcParams['legend.fontsize'] = 'small'
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['legend.borderpad'] = 0


import pandas as pd
import xlrd
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
sns.set_palette(custom_colors)
# sns.palplot(sns.color_palette("bright", 8))
# custom_colors = custom_colors[:1]+custom_colors[2:5]+custom_colors[6:7]+custom_colors[-1:]
# custom_colors = ["DarkBlue", "Crimson", "DarkGreen", "DarkOrchid", "Orange", "DeepPink"]
# You can use any HTML colors in the above line. For HTML colors, see this: https://www.w3schools.com/colors/colors_names.asp


# Setting current directory and target directory
user_folder = "F:\\Google Drive\\Research\\Projects"
# user_folder = "C:\\Users\\aravi\\Google Drive\\Research\\Projects"

##### Selecting the folder ####################################################
#dir_path = "\\Pd Interlayer contacts\\Cascade\\2019-10-11-I31\\"; Vg_idvd_plt = 60
#dir_path = "\\Pd Interlayer contacts\\Janis\\2019-11-07-Au\\"
#dir_path = "\\Pd Interlayer contacts\\Janis\\2019-11-21-I71\\"; Vg_idvd_plt = 20 
# dir_path = "\\InGaAs contacts\\Cascade\\2019-10-29-Ge-MoS2"; Vg_idvd_plt = 50
# dir_path = "\\InGaAs contacts\\Janis\\2019-10-31-Ge-MoS2"; Vg_idvd_plt = 50
# dir_path = "\\InGaAs contacts\\Janis\\2019-11-10-Ge-MoS2"; Vg_idvd_plt = 50
#dir_path = "\\Pd Interlayer contacts\\Janis\\2020-01-21-IA1\\"
# dir_path = "\\Pd Interlayer contacts\\Janis\\2020-01-31-IA3\\"; Vg_idvd_plt = 20
#dir_path = "\\PSG Doping\\Sentaurus\\P12_01 D2-4\\"
# dir_path = "\\InGaAs contacts\\Janis\\2020-01-29-ME2\\"; Vg_idvd_plt = 25
# dir_path = "\\InGaAs contacts\\Janis\\2020-02-18-ME5\\"; Vg_idvd_plt = 5
#dir_path = "\\Pd Interlayer contacts\\Janis\\2020-02-12-IA4\\"
# dir_path = "\\Pd Interlayer contacts\\Janis\\2020-02-27-IA3-ALOX\\"; Vg_idvd_plt = 30
#dir_path = "\\InGaAs contacts\\Cascade\\2020-03-03-MF1"; Vg_idvd_plt = 70
# dir_path = "\\Pd Interlayer contacts\\Janis\\2020-03-02-IB2\\"; Vg_idvd_plt = 20
# dir_path = "\\Pd Interlayer contacts\\Janis\\2020-08-30-IB13\\"; Vg_idvd_plt = 30
# dir_path = "\\Pd Interlayer contacts\\Janis\\2020-08-30-IB14\\"; Vg_idvd_plt = 30
# dir_path = "\\Pd Interlayer contacts\\Janis\\2020-09-02-IB13-ALOX\\"; Vg_idvd_plt = 30
# dir_path = "\\Pd Interlayer contacts\\Janis\\2020-09-02-IB14-ALOX\\"; Vg_idvd_plt = 30
# dir_path = "\\Pd Interlayer contacts\\Semi-Auto Cascade\\2020-09-01-IB13\\"; Vg_idvd_plt = 30
# dir_path = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-07-11-MF1"; Vg_idvd_plt = 40; Vds_low = 5;
# dir_path = "\\Pd Interlayer contacts\\Semi-Auto Cascade\\2020-09-01-IB13\\"; Vg_idvd_plt = 30
dir_path = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-19-IB21A"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
##############################################################################

target_dir = user_folder + dir_path
os.chdir(target_dir)  # Move to target directory

# Reading parameter files
# Reading channel parameters - number of channels, channel names, etc.
channel_params = Path('channel_parameters.xlsx')  # Channel parameters file
channel_params = xlrd.open_workbook(channel_params, logfile=open(os.devnull, 'w'))
channel_params = np.asarray(pd.read_excel(channel_params, engine='xlrd'))  # Reading Id-Vg data
channel_length = channel_params[:, 0]  # Channel lengths (string)
channel_L = channel_params[:, 1]  # Channel lengths (um)
channel_index = np.arange(channel_L.shape[0])  # Channel lengths index

# Constants
Cox_30nm = 116  # 30nm SiO2 capacitance (Kirby ACS Nano) (nF/cm2)

# Reading device parameters - file paths, column numbers, polarity of sweep, etc.
device_params = Path('device_parameters_auto.xlsx')  # Device parameters file - what sweeps, column indices, etc.
device_params = xlrd.open_workbook(device_params, logfile=open(os.devnull, 'w'))
device_params = np.asarray(pd.read_excel(device_params, engine='xlrd'))  # Reading Id-Vg data
device = device_params[:, 2]  # List of devices
file_prefix = device_params[:, 3]  # List of file prefixes
isfolder_organized = device_params[:, 4]  # Whether the device data is organized into folders
device_index = np.arange(device.shape[0])  # Device index

# Main section of code
os.chdir(target_dir)
tlm_set = [None] * device_index.shape[0]  # Creating None array to hold TLM objects

for i in device_index:
    if pd.isnull(file_prefix[i]):  # Empty file prefix is read as NaN from excel
        file_prefix[i] = ""
    if isfolder_organized[i]:  # If organized into folders, cd into folder
        os.chdir(device[i])
    inp_count = 0
    col_index = device_params[np.ix_([i], [5, 6, 7, 8, 9, 10, 19])]  # Column indices for various values
    col_index = col_index[0]  # Converting 2-D array to 1-D
    # 1,2,3 - Ids, Vgs, Vds for Id-Vg
    # 4,5,6 - Ids, Vgs, Vds for Id-Vd
    num_vd_idvg = device_params[i, 16]  # Number of Vd steps in Id-Vg sweep
    num_var_idvg = device_params[i, 17]  # Number of column variables in a single Id-Vg sweep
    num_vg_idvd = device_params[i, 11]  # Number of Vg points in Id-Vd sweep
    isbipolar_idvd = device_params[i, 12]  # Whether Id-Vd sweep is bipolar
    num_var_idvd = device_params[i, 13]  # Number of variables in a single Id-Vd measurement
    
    W = device_params[i, 14]  # Channel width (um)
    tox = device_params[i, 15]  # Thickness of the gate oxide (nm)
    Cox = Cox_30nm * 30 / tox  # 100nm SiO2 capacitance (nF/cm2)
    input_data = [None] * channel_length.shape[0]  # Creating dummy array for input
    
    channel_select = list(device_params[i, 18]) #Encoding to include/exclude channels
    
    for j in channel_index:        
        # filename_g = file_prefix[i] + channel_length[j] + "-G.csv"  # Filename of id-vg file
        # filename_d = file_prefix[i] + channel_length[j] + "-D.csv"  # Filename of id-vd file
        # filename_g = file_prefix[i] + channel_length[j] + "-G-AlOx.xls"  # Filename of id-vg file
        # filename_d = file_prefix[i] + channel_length[j] + "-D-AlOx.xls"  # Filename of id-vd file
        filename_g = file_prefix[i] + channel_length[j] + "-G-Vds0.1.csv"  # Filename of id-vg file
        filename_d = file_prefix[i] + channel_length[j] + "-D-Vds0.1.csv"  # Filename of id-vd file
        # Reading Id-Vg files
        my_file_g = Path(filename_g)  # Creating file object        
        my_file_d = Path(filename_d)  # Creating file object        
        
        if my_file_g.exists() and my_file_d.exists() and int(channel_select[j+1]):  # Checking if IdVg and IdVd files exist            
            print(my_file_g)
            wb_idvg = xlrd.open_workbook(filename_g, logfile=open(os.devnull, 'w'))
            idvg_data = np.asarray(pd.read_excel(wb_idvg, engine='xlrd'))  # Reading Id-Vg data            
            wb_idvd = xlrd.open_workbook(filename_d, logfile=open(os.devnull, 'w'))
            idvd_data = np.asarray(pd.read_excel(wb_idvd, engine='xlrd'))  # Reading Id-Vd data            
            input_data[inp_count] = InputData(idvg_data, channel_length, channel_L, j, idvd_data)  # Assigning the read data to InputData object
            inp_count += 1
    
    if inp_count:
        params = Parameters(col_index, num_var_idvg, num_vd_idvg, W, Cox, num_var_idvd, num_vg_idvd, isbipolar_idvd)
        
        tlm_set[i] = TLM(input_data[:inp_count], inp_count, params)  # Creating TLM object
        # tlm_set[i].tlm_calc(vd_step = 'low')  # Performing TLM calculations    
        tlm_set[i].tlm_calc(params, vd_step = 'hi')  # Performing TLM calculations    
        
        if tlm_set[i].count > 1:
            mpl.rcParams['figure.figsize'] = 12, 6
            fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3) #More than 1 channel works    
            # mpl.rcParams['figure.figsize'] = 16, 6
            # fig, ((ax1, ax2, ax3, ax7), (ax4, ax5, ax6, ax8)) = plt.subplots(2, 4) #More than 1 channel works    
        else:        
            mpl.rcParams['figure.figsize'] = 12, 3
            fig, ((ax1, ax2, ax3)) = plt.subplots(1, 3) #Only 1 channel works    
            # fig, ((ax1, ax2, ax3)) = plt.subplots(1, 3) #Only 1 channel works    
        
        ########### Plotting ############################################
        tlm_set[i].plot_IdVg(ax1, fig, lg = 'log', vd_step = 'low')  # Plotting Id-Vg for all channel lengths in log scale
        tlm_set[i].plot_IdVg(ax2, fig, lg = 'log', vd_step = 'hi')  # Plotting Id-Vg for all channel lengths in log scale
        tlm_set[i].plot_IdVg_multVds(ax3, fig, params, 'log', channel = 0)  # Plotting Id-Vg for 100nm in log scale for smallest and largest vds
        # tlm_set[i].plot_IdVd(ax2, fig, Vg_idvd_plt)  # Plotting Id-Vd for all channel lengths @ Vgs=50V
        # tlm_set[i].idvd[0].plot_IdVd(ax3, fig)  # Plotting Id-Vd for 1st channel length
        
        if tlm_set[i].count > 1:
            tlm_set[i].plot_TLM_fit(ax4, fig)  # Plotting TLM fit for the last Vov
            tlm_set[i].plot_Rc(ax5, fig, flag = 'n2D')  # Plotting Rc vs n2D
            tlm_set[i].plot_mu_TLM(ax6, fig, flag = 'n2D')  # Plotting mu_eff vs n2D
            # tlm_set[i].plot_Rsh(ax7, fig, 'n2D')  # Plotting Rsh vs n2D
            rsh = tlm_set[i].Rsh[-1]
            rc = tlm_set[i].Rc[-1]
            
        fig.suptitle(device[i], fontsize=14)
        plt.show()
        fig.savefig("summary-" + str(i) + ".svg", transparent=True,\
                                        bbox_inches='tight', pad_inches=0.1)    
    if isfolder_organized[i]:
        os.chdir(target_dir)