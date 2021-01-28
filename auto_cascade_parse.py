"""
Created on Tue March 5 03:14:00 2020
Updated on Sunday 12/06/2020
@author: Aravindh Kumar
"""
""" 
Python3 code to:
1) Rename auto-cascade files in a directory or folder 
to make it easier to parse for analysis

2) Generate yield summary plots

3) Generate device parameter and channel parameter files to analyze auto-cascade data
"""

# Importing modules
import os 
from shutil import copyfile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

plt.close('all')

#Changing plot properties
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
mpl.rcParams['figure.figsize'] = 6, 6  # Figure size
mpl.rcParams['figure.constrained_layout.h_pad'] = 0.04167  # Padding around axes objects
mpl.rcParams['figure.constrained_layout.w_pad'] = 0.04167
mpl.rcParams['legend.fontsize'] = 'small'
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['legend.borderpad'] = 0
mpl.rcParams['image.aspect'] = 0.5


def heat_map(ax, Z):
    xlabels = [str(i) for i in range(1,Z.shape[1]+1)]
    ylabels = [str(i) for i in range(1,Z.shape[0]+1)]        
    #Major ticks
    ax.set_xticks(np.arange(len(xlabels)))
    ax.set_yticks(np.arange(len(ylabels)))    
    #Major tick labels
    ax.set_xticklabels(xlabels)
    ax.set_yticklabels(ylabels)    
    #Minor ticks
    ax.set_xticks(np.arange(-0.5, Z.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, Z.shape[0], 1), minor=True)    
    ax.imshow(Z, cmap='binary_r')   
    ax.grid(which='minor', color='b', linestyle='-', linewidth=2)

# Function to rename multiple files 
def rename(user_folder, target_dir, num_rows = 15, num_cols = 9):
    """ 
    Case: 1
    Input - IDVG_Die_X5_Y-1_Dev3_Rep1.csv
    Output - 1-5-CH3-G.csv

    Case: 2
    Input - IDVG_Die_X5_Y-15_Dev3_Rep1_Vds0.1.csv / IDVG_Die_X5_Y-15_Dev3_Rep1_Vds1.csv
    Output - 15-5-CH3-G-Vds0.1.csv / 15-5-CH3-G-Vds1.csv
    """
    
    source_dir = user_folder + target_dir + "Source Data\\"
    target_dir = user_folder + target_dir
    
    gate_leak = np.zeros((num_rows,num_cols))
    working_device = np.zeros((num_rows,num_cols))
    open_device = np.zeros((num_rows,num_cols))
    short_device = np.zeros((num_rows,num_cols))
    low_curr = np.zeros((num_rows,num_cols))    

    for filename in os.listdir(source_dir):     
    # for filename in os.listdir(target_dir):   
        index_working = filename.find('IDVG')
        index_open = filename.find('Open')
        index_gateleak = filename.find('GateLeakage')
        index_short = filename.find('Short')
        index_lowcurr = filename.find('LowCurr')
        index1 = filename.find('_Die_X')
        index2 = filename.find('_Y-')
        index3 = filename.find('_Dev')
        index4 = filename.find('.csv')
        index5 = filename.find('_Vds')
#        print("Index1 ="+str(index1)+" Index2 ="+str(index2)+" Index1 ="+str(index3))        
        if index1 != -1 and index2 != -1 and index3 != -1 and index4 != -1:
            row_index = int(filename[index2+3:index3])
            ch_index = int(filename[index3+4])
            col_index = int(filename[index1+6])
            if index_working != -1:            
                working_device[row_index-1, col_index-1] = 1
    #            print("Row index ="+str(row_index)+" Col index ="+str(col_index)+" Ch index ="+str(ch_index))
                if index5 != -1:
                    new_filename = str(row_index)+"-"+str(col_index)+"-CH"+str(ch_index)+"-G-"+filename[index5+1:]
                else:
                    new_filename = str(row_index)+"-"+str(col_index)+"-CH"+str(ch_index)+"-G"+filename[index4:]
    #            print(new_filename)        
                src = source_dir +'\\'+ filename
                # src = target_dir +'\\'+ filename    
                dst = target_dir +'\\'+ new_filename    
                # rename() function will 
                # rename all the files 
                copyfile(src, dst)
            elif index_open != -1:
                open_device[row_index-1, col_index-1] += 1
            elif index_gateleak != -1:
                gate_leak[row_index-1, col_index-1] = 1
            elif index_short != -1:
                short_device[row_index-1, col_index-1] = 1
            elif index_lowcurr != -1:
                low_curr[row_index-1, col_index-1] = 1
    open_device = (open_device >= 4).astype(int)
    
    return gate_leak, working_device, open_device, short_device, low_curr

def plot_heatmaps(gate_leak, working_device, open_device, short_device, low_curr):
    """ Plot the heatmaps to show the map of leaking, working, open, short and low-current dies"""

    fig1, ax1 = plt.subplots(figsize=(6, 6))    
    heat_map(ax1, gate_leak)
    fig1.suptitle("Gate Leakage")
    fig1.savefig('gateleak.svg',transparent=True, bbox_inches='tight', pad_inches=0.1)
    
    fig2, ax2 = plt.subplots()        
    heat_map(ax2, working_device)
    fig2.suptitle("Working Devices")
    fig2.savefig('working.svg',transparent=True, bbox_inches='tight', pad_inches=0.1)
    
    fig3, ax3 = plt.subplots()        
    heat_map(ax3, open_device)
    fig3.suptitle("Open Devices")
    fig3.savefig('open.svg',transparent=True, bbox_inches='tight', pad_inches=0.1)
    
    fig4, ax4 = plt.subplots()        
    heat_map(ax4, short_device)
    fig4.suptitle("Shorted Devices")
    fig4.savefig('short.svg',transparent=True, bbox_inches='tight', pad_inches=0.1)
    
    fig5, ax5 = plt.subplots()        
    heat_map(ax5, low_curr)
    fig5.suptitle("Low-current Devices")
    fig5.savefig('low_curr.svg',transparent=True, bbox_inches='tight', pad_inches=0.1)
    
    return fig1, fig2, fig3, fig4, fig5


def create_device_parameters(user_folder, target_dir, W, tox, num_var_IdVg, num_Vd_IdVg, num_rows = 15, num_cols = 9):
    """
    Generate device parameter file for auto-cascade data - device_parameters_auto.xlsx
    """
    num_device = 0
    device_number = []
    device = []
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

    path = user_folder + target_dir

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
    # channel_lengths = [0.1, 4, 0.3, 5, 0.7, 1]
    # channel_labels = ["CH1", "CH2", "CH3", "CH4", "CH5"]
    # channel_lengths = [4, 5, 6, 7, 10]
    channel_data = {"Channel Lengths":channel_labels, "Channel L":channel_lengths}
    df = pd.DataFrame(channel_data)
    df.to_excel(path+"channel_parameters.xlsx",index=False)

user_folder = "F:\\Google Drive\\Research\\Projects"
# user_folder = "C:\\Users\\aravi\\Google Drive\\Research\\Projects"

# target_dir = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-03-04-MF1\\"; W = 10; tox = 100; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-07-11-MF1\\"; W = 10; tox = 100; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-02-MH3\\"; W = 1; tox = 95;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-13-MJ2\\"; W = 1; tox = 95;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-13-MJ3\\"; W = 1; tox = 95;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-23-MJ1\\"; W = 1; tox = 95;
# target_dir = "\\InGaAs Contacts\\Semi-Auto Cascade\\2020-10-29-MJ5\\"; W = 1; tox = 95;
# target_dir = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-07-11-MF1\\"; W = 10; tox = 100; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-12-11-MJ6\\"; W = 1; tox = 95; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\InGaAs contacts\\Semi-Auto Cascade\\2020-12-11-MJ6-Large\\"; W = 1; tox = 95; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-07-22-IB3\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-07-22-IB4\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-07-22-IB5\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-07-27-IB7\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-07-27-IB8\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-07-27-IB9\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-01-IB13\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-01-IB14\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-09-06-IB13-source-ground\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-01-IB13-N2-anneal\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-03-IB13-N2-anneal\\"; W = 2; tox = 88;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-10-20-IB18-N2-purge\\"; W = 1; tox = 95;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-19-IB21A\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-19-IB21B\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-24-IB21A-capped\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-05-IB21A-capped-longer-sweep\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-22-IB21A-IdVds\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-07-IB21B-capped\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-29-IB22A\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-11-29-IB22B\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-03-IB22A-annealed\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-07-IB22A-capped\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-07-IB22B-capped\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;
# target_dir = "\\Pd Interlayer Contacts\\Semi-Auto Cascade\\2020-12-12-IB23\\"; W = 2; tox = 50; num_var_IdVg = 4; num_Vd_IdVg = 1;


num_rows = 15
num_cols = 9

gate_leak, working_device, open_device, short_device, low_curr = rename(user_folder, target_dir, num_rows, num_cols)
# plot_heatmaps(gate_leak, working_device, open_device, short_device, low_curr)
create_device_parameters(user_folder, target_dir, W, tox, num_var_IdVg, num_Vd_IdVg, num_rows, num_cols)

plt.close('all')