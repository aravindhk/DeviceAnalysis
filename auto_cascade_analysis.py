"""
Created on Thursday September 17 20:15:00 2020
Updated on Sunday 12/06/2020
@author: Aravindh Kumar
"""

"""
Python3 code to:
Analyse auto-cascade data
Currently only Id-Vg analysis is supported
"""

import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import xlrd
import numpy as np
import seaborn as sns
from pathlib import Path
from device_analysis_classes import *
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings("ignore")

# Changing plot properties
mpl = plot_properties(mpl)

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

interp = 0 #Specify whether Id-Vg interpolation required
gateLeak = 0 #Number of gate leakage sites
calc_flag = 1 #Specify whether TLM calculations should be carried out
decimals_Vg = 0 #Number of decimals in Vgs after rounding
# Constants
Cox_30nm = 116 #30nm SiO2 capacitance (Kirby ACS Nano) (nF/cm2)

flags = Flags(num_rows = 15, num_cols = 9)    

# Setting current directory and target directory
user_folder = "F:\\Google Drive\\Research\\Projects"
# user_folder = "C:\\Users\\aravi\\Google Drive\\Research\\Projects"
##### Selecting the folder ####################################################
# target_dir = "\\WS2 Contacts and Doping\\Janis\\2019-03-09_TLMs_on_WS2\\Before AlOx"; Vg_idvd_plt = 80; Vds_low = 5.0; isbipolar_idvg = 0; decimals_Vg = 1;
# target_dir = "\\WS2 Contacts and Doping\\Janis\\2019-03-09_TLMs_on_WS2\\After AlOx"; Vg_idvd_plt = 80; Vds_low = 5.0; isbipolar_idvg = 0; decimals_Vg = 1;
# target_dir = "\\WS2 Contacts and Doping\\Janis\\2019-03-19-TLM-WS2-RUN2"; Vg_idvd_plt = 80; Vds_low = 1.0; isbipolar_idvg = 0; decimals_Vg = 1;
# target_dir = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-03-04-MF1\\"; Vg_idvd_plt = 40; Vds_low = 5; isbipolar_idvg = 1;
# target_dir = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-07-11-MF1\\"; Vg_idvd_plt = 40; Vds_low = 5; isbipolar_idvg = 1; decimals_Vg = 2;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-02-MH3\\"; Vg_idvd_plt = 30; Vds_low = 0.5; isbipolar_idvg = 1;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-13-MJ2\\"; Vg_idvd_plt = 30; Vds_low = 0.5; isbipolar_idvg = 1;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-13-MJ3\\"; Vg_idvd_plt = 30; Vds_low = 0.5; isbipolar_idvg = 1;
# target_dir = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-12-11-MJ6\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-12-11-MJ6-Large\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-01-31-IA3"; Vg_idvd_plt = 20; Vds_low = 0.1; isbipolar_idvg = 0;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-01-IB13"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 0;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-01-IB14"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 0;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-06-IB13-source-ground"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-01-IB13-N2-anneal"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-03-IB13-N2-anneal"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-20-IB18-N2-purge\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-19-IB21A\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-19-IB21B\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-24-IB21A-capped\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-05-IB21A-capped-longer-sweep\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1; interp = 1;
target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-22-IB21A-IdVds\\"; Vg_idvd_plt = 30; Vds_low = 1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-07-IB21B-capped\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-29-IB22A\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-29-IB22B\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-03-IB22A-annealed\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-07-IB22A-capped\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-07-IB22B-capped\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-12-IB23\\"; Vg_idvd_plt = 30; Vds_low = 0.1; isbipolar_idvg = 1;
###############################################################################
target_dir = user_folder + target_dir
os.chdir(target_dir) #Move to target directory

#Reading channel parameters and Device parameters
channel_params, device_params, channel_length, channel_L, channel_label, channel_index, chip,\
    device, file_prefix, isfolder_organized, W, tox, num_channels = read_params(target_dir)
Cox = Cox_30nm*30/tox #100nm SiO2 capacitance (nF/cm2)
device_index = np.arange(device.shape[0]) #Device index

# Main section of code
os.chdir(target_dir)
tlm_set = [None]*device_index.shape[0] #Creating None array to hold TLM objects
pp = PdfPages('foo.pdf')
device_count = 0

for i in device_index:        
    if pd.isnull(file_prefix[i]): #Empty file prefix is read as NaN from excel
        file_prefix[i] = ""
    if isfolder_organized[i]: #If organized into folders, cd into folder
        os.chdir(chip[i])
            
    col_index = device_params[np.ix_([i],[5,6,7,8,9,10,19])] #Column indices for various values
    col_index = col_index[0] #Converting 2-D array to 1-D
    #1,2,3 - Ids, Vgs, Vds for Id-Vg
    #4,5,6 - Ids, Vgs, Vds for Id-Vd    
    num_vd_idvg = device_params[i, 16]  # Number of Vd steps in Id-Vg sweep
    num_var_idvg = device_params[i, 17]  # Number of column variables in a single Id-Vg sweep
    # isbipolar_idvd = device_params[i, 12]  # Whether Id-Vd sweep is bipolar
    # num_vg_idvd = device_params[i, 11]  # Number of Vg steps in Id-Vd sweep
    # num_var_idvd = device_params[i, 13]  # Number of column variables in a single Id-Vd sweep

    W = device_params[i, 14]  # Channel width (um)
    tox = device_params[i, 15]  # Thickness of the gate oxide (nm)
    Cox = Cox_30nm * 30 / tox  # 100nm SiO2 capacitance (nF/cm2)
    input_data = [None] * channel_length.shape[0]  # Creating dummy array for input    
    channel_select = list(device_params[i, 18]) #Encoding to include/exclude channels
    channel_count = 0    
    
    for j in channel_index:    
        filename_g = file_prefix[i] + channel_length[j] + "-G-Vds" + str(Vds_low) #Filename of id-vg files        
        # filename_g = file_prefix[i] + channel_length[j] + "-G" #Filename of id-vg files        
        # filename_g = file_prefix[i] + channel_length[j] + "-G-1"  # Filename of id-vg file
        # filename_g = file_prefix[i] + channel_length[j]  # Filename of id-vg file        
        # filename_d = file_prefix[i] + channel_length[j] + "-D" #Filename of id-vd file
        #Reading Id-Vg files
        file_g_xls = Path(filename_g + ".xls") #Creating file object
        file_g_csv = Path(filename_g + ".csv") #Creating file object        
        # my_file_d = Path(filename_d) #Creating file object
        exist_flag = 0        
        if file_g_csv.exists() and int(channel_select[j+1]): #Checking if IdVg and IdVd files exist                    
            idvg_data = csv_read(filename_g)
            exist_flag = 1
        elif file_g_xls.exists() and int(channel_select[j+1]): #Checking if IdVg and IdVd files exist
            idvg_data = xls_read(filename_g)
            exist_flag = 1
        if exist_flag == 1:                        
            input_data[channel_count] = InputData(idvg_data, channel_length, channel_L, j) #Assigning the read data to InputData object
            channel_count += 1                 
            
    if channel_count:           
        print(device[i])
        row, column = device[i].split('-')        
        params = Parameters(col_index, num_var_idvg, num_vd_idvg, W, Cox,\
                            isbipolar_idvg = isbipolar_idvg, Vds_param = Vds_low,\
                                decimals_Vg = decimals_Vg, interp = interp, row = row, column = column)
        # tlm_set[device_count] = TLM(input_data[:channel_count],channel_count,\
        #                             params, column = int(column), row = int(row)) #Creating TLM object
        tlm_set[device_count] = TLM(input_data[:channel_count],channel_count,\
                                    params) #Creating TLM object
        if 0:#tlm_set[device_count].gateLeak:
            flags.gate_leak[int(row)-1, int(column)-1] = 1
            gateLeak = gateLeak + 1
        if calc_flag == 1:
            tlm_set[device_count].tlm_calc() #Performing TLM calculations                    
        if tlm_set[device_count].R_squared>=0.0 and tlm_set[device_count].count>2 and tlm_set[device_count].goodTLM:    
            print("Good TLM " + str(tlm_set[device_count].R_squared))
            fig_summary, ax_summary = plt.subplots(2, 3) #More than 1 channel works
            ax_summary = ax_summary.reshape(-1)
            tlm_set[device_count].summary_plot(ax_summary, fig_summary, device[i])
            pp.savefig(fig_summary)
            if tlm_set[device_count].R_squared>=0.9:
                fig_summary.savefig('summary'+tlm_set[device_count].position+'.svg',\
                                    transparent=True, bbox_inches='tight', pad_inches=0.1)
            plt.close(fig_summary)
            if isfolder_organized[i]:
                os.chdir(target_dir)                    
        device_count = device_count + 1 #To keep count of analyzed TLMs
    
    if isfolder_organized[i]:
        os.chdir(target_dir)

fig, ax = plt.subplots(2,3) #To plot Id-Vg by channel
ax = ax.reshape(-1)
plot_by_Lch(fig, ax, device_count, tlm_set, channel_L, channel_label, num_channels)
fig.savefig('summary.svg',transparent=True, bbox_inches='tight', pad_inches=0.1)
pp.savefig(fig)

fig_gate_leak_map = plot_heatmaps(gate_leak = flags.gate_leak)
# fig_gate_leak_map.savefig('gateleak_map.svg',transparent=True, bbox_inches='tight', pad_inches=0.1, dpi=300)
pp.savefig(fig_gate_leak_map)

fig_gate, ax_gate = plt.subplots(figsize=(4, 3))
plot_gate_leakage(fig_gate, ax_gate, device_count, tlm_set, channel_L, channel_length)
fig_gate.savefig('gateleak.svg',transparent=True, bbox_inches='tight', pad_inches=0.1, dpi=300)
pp.savefig(fig_gate)

fig_rtot, ax_rtot = plt.subplots(figsize=(4, 3)) #To plot Id-Vg by channel
plot_all_Rtot(fig_rtot, ax_rtot, device_count, tlm_set, n2D = 1e13)
fig_rtot.savefig('rtot.svg',transparent=True, bbox_inches='tight', pad_inches=0.1, dpi=300)
pp.savefig(fig_rtot)

fig_rc_cdf, ax_rc_cdf = plt.subplots(figsize=(4, 3)) #To plot Id-Vg by channel
plot_Rc_cdf(fig_rc_cdf, ax_rc_cdf, device_count, tlm_set)
fig_rc_cdf.savefig('rc_cdf.svg',transparent=True, bbox_inches='tight', pad_inches=0.1, dpi=300)
pp.savefig(fig_rc_cdf)

pp.close()
plt.close('all')          