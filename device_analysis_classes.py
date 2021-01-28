"""
Created on Thursday September 17 18:30:00 2020
Updated on Sunday 12/06/2020
@author: Aravindh Kumar
"""

"""
Python3 code:
Class definitions to analyse Cascade/Janis data
"""


import os
import matplotlib.pyplot as plt
from scipy.signal import find_peaks  # For finding local maxima in gm/dgm

import pandas as pd
import xlrd
import numpy as np
import seaborn as sns
import statsmodels.api as sm
from pathlib import Path
from matplotlib import ticker


# Constants
q = 1.6e-19  # Charge of an electron (C)
Id_thresh = 2e-8  # For constant Vt calculation (A)
I_noise = 1e-12  # Noise floor for Id (A)
n2D_LL = 1e10  # n2D lower limit for plotting (cm-2)
compliance_threshold = 1e2 # Threshold to determine whether device hit compliance

# Plot labels
# vd_label = r'Drain Voltage,  $V_{DS}$  [V]'
# vg_label = r'Gate Voltage,  $V_{GS}$  [V]'
# vov_label = r'Overdrive Voltage,  $V_{ov}$  [V]'
# n2d_label = r'Carrier Concentration,  $n_{2D}$  [$10^{12}$ $cm^{-2}$]'
# id_label = r'Drain Current,  $I_D$  [$\mu A$/$\mu m$]'
# ig_label = r'Gate Current,  $I_G$  [$A$]'
# rc_label = r'Contact Resistance,  $R_{C}$  [k$\Omega\cdot\mu m$]'
# mu_eff_label = r'Effective Mobility,  $\mu_{eff}$  [$cm^2$$V^{-1}$$s^{-1}$]'
# mu_fe_label = r'Field-effect Mobility, ${\mu}_{FE}$ [$cm^2$$V^{-1}$$s^{-1}$]'
# gm_label = r'Transconductance $g_{m}$ [$S$/$\mu$$m$]'
# rsh_label = r'Sheet Resistance, $R_{sh}$ [k$\Omega$\cdot\boxdot]'
# rtot_label = r'Total Resistance, $R_{tot}$ [k$\Omega\cdot\mu $m]'
# lch_label = r'Channel Length,  $L_{ch}$  [$\mu$m]'

vd_label = r'$V_{DS}$  [V]'
vg_label = r'$V_{GS}$  [V]'
vov_label = r'$V_{ov}$  [V]'
n2d_label = r'$n_{2D}$  [$10^{12}$ $cm^{-2}$]'
id_label = r'$I_D$  [$\mu A$/$\mu m$]'
ig_label = r'$I_G$  [$A$]'
rc_label = r'$R_{C}$  [k$\Omega\cdot\mu m$]'
rc_cdf_label = r'$R_{C}$  [$\Omega\cdot\mu m$]'
mu_eff_label = r'$\mu_{eff}$  [$cm^2$$V^{-1}$$s^{-1}$]'
mu_fe_label = r'${\mu}_{FE}$ [$cm^2$$V^{-1}$$s^{-1}$]'
gm_label = r'$g_{m}$ [$S$/$\mu$$m$]'
rsh_label = r'$R_{sh}$ [k$\Omega$\cdot\boxdot]'
rtot_label = r'$R_{tot}$ [k$\Omega\cdot\mu $m]'
lch_label = r'$L_{ch}$  [$\mu$m]'

annotate_x = 0.05
annotate_y1 = 0.85
annotate_y2 = 0.75
annotate_y22 = 0.725


# Helper functions
def csv_read(filename_g):
    """ Used to read input data in csv format"""    

    filename_g = filename_g + ".csv"
    file_g = Path(filename_g)        
    idvg_data = np.asarray(pd.read_csv(file_g))
    return idvg_data


def xls_read(filename_g):
    """ Used to read input data in xls format"""   

    filename_g = filename_g + ".xls"
    idvg_data = np.asarray(pd.read_excel(xlrd.open_workbook(filename_g, logfile=open(os.devnull, 'w')), engine='xlrd'))  # Reading Id-Vg data
    return idvg_data


def plot_properties(mpl):
    """ Sets the matplotlib rcParams """

    mpl.rcParams['lines.linewidth'] = 1
    mpl.rcParams['axes.titlesize'] = 10
    mpl.rcParams['xtick.labelsize'] = 9
    mpl.rcParams['ytick.labelsize'] = 9
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['font.size'] = 10
    mpl.rcParams['lines.markersize'] = 3
    mpl.rcParams['figure.autolayout'] = True
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
    mpl.rcParams["errorbar.capsize"] = 2
    return mpl


def read_params(target_dir, channel_params_file = "channel_parameters.xlsx", device_params_file = "device_parameters_auto.xlsx"):
    """ Read the parameters from channel_params_file and device_params_file"""
    
    os.chdir(target_dir)

    #Reading channel parameters - how many channels, what channel names, etc.
    channel_params = Path(channel_params_file) #Channel parameters file
    channel_params = xlrd.open_workbook(channel_params,logfile=open(os.devnull,'w'))
    channel_params = np.asarray(pd.read_excel(channel_params,engine='xlrd')) #Reading Id-Vg data
    
    #Reading device parameters
    device_params = Path(device_params_file) #Device parameters file - type of sweep, column indices, etc.
    device_params = xlrd.open_workbook(device_params,logfile=open(os.devnull,'w'))
    device_params = np.asarray(pd.read_excel(device_params,engine='xlrd')) #Reading Id-Vg data

    #Extracting the parameters    
    channel_length = channel_params[:,0] #Channel lengths (string)
    channel_L = channel_params[:,1] #Channel lengths (um)
    channel_label = np.empty(channel_L.shape[0], dtype=object)        
    num_channels = channel_L.shape[0]
    for i in np.arange(channel_L.shape[0]):
        if channel_L[i] < 1:
            channel_label[i] = str(int(channel_L[i]*1000))+r' nm'
        else:
            channel_label[i] = str(channel_L[i])+r' $\mu$m'
    
    channel_index = np.arange(channel_L.shape[0]) #Channel lengths index
    chip = device_params[:,1] #List of chips (needed?)
    device = device_params[:,2] #List of devices
    file_prefix = device_params[:,3] #List of file prefixes
    isfolder_organized = device_params[:,4] #Whether the device data is organized into folders
    W = device_params[0,14] #Channel width (um)
    tox = device_params[0,15] #Thickness of the gate oxide (nm)

    return channel_params, device_params, channel_length, channel_L,\
        channel_label, channel_index, chip, device, file_prefix, isfolder_organized, W, tox, num_channels
    

def cdf(fig, ax, x, xerr = np.empty(shape=(0), dtype=object)):
    """CDF plot for Rc"""

    sort_index = np.argsort(x)
    x, y = sorted(x), np.arange(len(x)) / (len(x) - 1)
    
    if xerr.shape[0] == 0:
        ax.plot(x, y)
    else:
        xerr = xerr[sort_index]
        ax.errorbar(x,y,xerr = xerr, ecolor='#87b5ff')
    ax.grid()


def smooth(a, WSZ=7):
    """ a: NumPy 1-D array containing the data to be smoothed
    WSZ: smoothing window size needs, which must be odd number,
    as in the original MATLAB implementation"""

    out0 = np.convolve(a, np.ones(WSZ, dtype=int), 'valid') / WSZ
    r = np.arange(1, WSZ - 1, 2)
    start = np.cumsum(a[:WSZ - 1])[::2] / r
    stop = (np.cumsum(a[:-WSZ:-1])[::2] / r)[::-1]
    return np.concatenate((start, out0, stop))


def find_nearest(array, value):
    """Returns the element closest to value in the array"""

    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def set_x_ticks(ax, num_ticks = 5,round_factor = None):
    start, end = ax.get_xlim()
    plot_range = np.abs(end - start)
    if round_factor == None:
        if plot_range >= 350:
            round_factor = 100
        # elif plot_range >= 250:
        #     round_factor = 75
        elif plot_range >= 200:
            round_factor = 50
        elif plot_range >= 80:
            round_factor = 20
        elif plot_range >= 60:
            round_factor = 15
        elif plot_range >= 40:
            round_factor = 10
        elif plot_range >= 20:
            round_factor = 5
        else:
            round_factor = 1
    tick_frequency = round((plot_range/num_ticks)/round_factor)*round_factor    
    ax.xaxis.set_ticks(np.arange(round(start/tick_frequency)*tick_frequency,\
                 (round(end/tick_frequency)+1)*tick_frequency, tick_frequency))
          
        
def set_ticks(ax, num_ticks = 5,axis = 'y',round_factor = None):
    """Set tick frequency"""

    if axis == 'y':
        start, end = ax.get_ylim()
        plot_range = np.abs(end - start)
    else:
        start, end = ax.get_xlim()
        plot_range = np.abs(end - start)
        
    if round_factor == None:
        if plot_range >= 350:
            round_factor = 100
        # elif plot_range >= 250:
        #     round_factor = 75
        elif plot_range >= 200:
            round_factor = 50
        elif plot_range >= 80:
            round_factor = 20
        elif plot_range >= 40:
            round_factor = 10
        elif plot_range >= 20:
            round_factor = 5
        else:
            round_factor = 1
            
    if axis == 'y':        
        tick_frequency = round((plot_range/num_ticks)/round_factor)*round_factor
        ax.yaxis.set_ticks(np.arange(round(start/tick_frequency)*tick_frequency,\
                 (round(end/tick_frequency)+1)*tick_frequency, tick_frequency))
    else:
        tick_frequency = round((plot_range/num_ticks)/round_factor)*round_factor
        ax.xaxis.set_ticks(np.arange(round(start/tick_frequency)*tick_frequency,\
                 (round(end/tick_frequency)+1)*tick_frequency, tick_frequency))
          
            
def heat_map(ax, Z, colormap = 'binary_r'):
    xlabels = [str(i) for i in range(1,Z.shape[1]+1)]
    ylabels = [str(i) for i in range(1,Z.shape[0]+1)]        
    #Major ticks
    ax.set_xticks(np.arange(len(xlabels)))
    ax.set_yticks(np.arange(len(ylabels)))    
    #Major tick labels
    ax.set_xticklabels(xlabels)
    ax.set_yticklabels(ylabels)    
    #Minor ticks
    ax.set_xticks(np.arange(-.5, Z.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, Z.shape[0], 1), minor=True)    
    ax.imshow(Z, cmap=colormap,aspect='auto')   
    ax.grid(which='minor', color='w', linestyle='-', linewidth=1)
    
    
def plot_heatmaps(gate_leak=None, working_device=None, open_device=None,\
                                          short_device=None, low_curr=None):
    """ Plot the heatmaps to show the map of leaking, 
    working, open, short and low-current dies"""

    if gate_leak is not None:
        fig1, ax1 = plt.subplots(figsize=(4, 4))    
        heat_map(ax1, gate_leak,colormap = 'seismic')
        ax1.set_title("Gate Leakage")
        fig1.savefig('gateleak.svg',transparent=True, bbox_inches='tight', \
                     dpi=300, pad_inches=0.1)
    
    if working_device is not None:
        fig2, ax2 = plt.subplots()        
        heat_map(ax2, working_device)
        ax2.set_title("Working Devices")
        fig2.savefig('working.svg',transparent=True, bbox_inches='tight', \
                     dpi=300, pad_inches=0.1)        
    
    if open_device is not None:
        fig3, ax3 = plt.subplots()        
        heat_map(ax3, open_device)
        ax3.set_title("Open Devices")
        fig3.savefig('open.svg',transparent=True, bbox_inches='tight', \
                     dpi=300, pad_inches=0.1)
    
    if short_device is not None:
        fig4, ax4 = plt.subplots()        
        heat_map(ax4, short_device)
        ax4.set_title("Shorted Devices")
        fig4.savefig('short.svg',transparent=True, bbox_inches='tight', \
                     dpi=300, pad_inches=0.1)
    
    if low_curr is not None:
        fig5, ax5 = plt.subplots()        
        heat_map(ax5, low_curr)
        ax5.set_title("Low-current Devices")
        fig5.savefig('low_curr.svg',transparent=True, bbox_inches='tight', \
                     dpi=300, pad_inches=0.1)


def plot_by_Lch(fig, ax, device_count, tlm_set, channel_L, channel_label, num_channels):
    """Plot the Id-Vg for all TLMs with each channel on a separate subplot"""
    
    for l in np.arange(num_channels):
        ax[l].annotate(r'$L_{CH}$ = ' + channel_label[l], xy=(annotate_x,\
                                       annotate_y1), xycoords='axes fraction')        
        ax[l].set_ylim(1e-6, 1e3)  # Sets lower y-limit        
    for i in np.arange(device_count):
        for k in np.arange(tlm_set[i].count):
            ax_index = np.where(channel_L == tlm_set[i].L[k])[0]
            ax_index = ax_index[0]
            if tlm_set[i].idvg[k,0].gateLeak == 0:
                tlm_set[i].plot_IdVg_multVds(ax[ax_index], fig, 'log', channel = k, anno = 'no')
                
                
def plot_cdf(fig, ax, device_count, tlm_set, get_function, n2D = 1e13):
    """Plots the value returned by get_function in a CDF plot at the specifed n2D for all TLMs"""

    X_array = np.zeros(device_count)
    del_X_array = np.zeros(device_count)
    R_squared_array = np.zeros(device_count)
    position_array = np.empty(shape=(device_count,), dtype='object')
    goodTLM_count = 0
    for i in np.arange(device_count):
        if tlm_set[i].gateLeak == 0 and tlm_set[i].goodTLM and tlm_set[i].R_squared >= 0.8:
            X, del_X, R_squared, goodX, xaxis_label = get_function(tlm_set[i], n2D = n2D)
            position = tlm_set[i].position
            if goodX:
                X_array[goodTLM_count] = X
                del_X_array[goodTLM_count] = del_X
                R_squared_array[goodTLM_count] = R_squared
                position_array[goodTLM_count] = position,
                goodTLM_count = goodTLM_count + 1
    rc_dataset = pd.DataFrame({'x': X_array[:goodTLM_count], \
        'del_x': del_X_array[:goodTLM_count], 'r_squared': R_squared_array[:goodTLM_count],\
            'position': position_array[:goodTLM_count]})
    rc_dataset.to_csv('rc.csv')
    print(goodTLM_count)
    cdf(fig, ax, X_array[:goodTLM_count], del_X_array[:goodTLM_count])
    ax.set(xlabel=xaxis_label, ylabel='CDF')


def plot_all_Rtot(fig, ax, device_count, tlm_set, n2D = 1e13):
    """Plots Rtot vs Lch for all TLMs at the specified n2D"""

    # ax.set_xlim(0.0,1.1)
    ax.set_ylim(0,50)    
    # ax.set_ylim(0,2e6)
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((0,0)) 
    ax.yaxis.set_major_formatter(formatter) 
    ax.set(xlabel=lch_label, ylabel=rtot_label)
    n2D_anno = n2D / 1e12
    if n2D_anno < 10:
        ax.annotate(r'at $n_{2D}$ = ' + str(np.round(n2D_anno, 1)) + r' $\times$ $10^{12} cm^{-2}$',
            xy=(annotate_x, annotate_y22), xycoords='axes fraction', fontsize=8.5)      
    else:
        ax.annotate(r'at $n_{2D}$ = ' + str(np.round(n2D_anno/10, 1)) + r' $\times$ $10^{13} cm^{-2}$',
            xy=(annotate_x, annotate_y22), xycoords='axes fraction', fontsize=8.5)                          
    # ax.set_yscale('log')
    # ax.set_ylim(1e3,1e8)    
    for i in np.arange(device_count):        
        if tlm_set[i].gateLeak == 0 and tlm_set[i].goodTLM:
            tlm_set[i].plot_Rtot(ax, fig, n2D = n2D, anno = 'no')            
    ax.grid()
            

def plot_gate_leakage(fig_gate, ax_gate, device_count, tlm_set, channel_L, channel_length):
    """Plot the gate for all TLMs for all channels on a single plot"""
    
    for i in np.arange(device_count):
        for k in np.arange(tlm_set[i].count):
            tlm_set[i].idvg[k,0].plot_Ig(ax_gate, fig_gate)  # Plotting Ig for 100nm in log scale for smallest Vds
    start, end = ax_gate.get_ylim()
    ax_gate.set_ylim(bottom=np.maximum(1e-12, start))  # Sets lower y-limit
    ax_gate.set_ylim(top=np.minimum(1e1, end))  # Sets upper y-limit        
    start, end = ax_gate.get_xlim()
    ax_gate.set_xlim(left=np.maximum(-30, start))  # Sets lower x-limit
    ax_gate.set_xlim(right=np.minimum(30, end))  # Sets upper x-limit        
    
        
class InputData:
    """Class for input data"""

    def __init__(self, idvg, channel_length, channel_L, index, idvd = None):
        self.id_vg = idvg  # Assigning read Id-Vg data
        self.channel_length = channel_length[index]  # Channel length (string)
        self.channel_L = channel_L[index]  # Channel length (value)
        self.id_vd = idvd  # Assigning read Id-Vd data


class Flags:
    """Class to store flags for gate leakage, working devices, etc"""
 
    def __init__(self, num_rows, num_cols):
        self.gate_leak = np.zeros((num_rows,num_cols))
        self.working_device = np.zeros((num_rows,num_cols))
        self.open_device = np.zeros((num_rows,num_cols))
        self.short_device = np.zeros((num_rows,num_cols))
        self.low_curr = np.zeros((num_rows,num_cols))  
    

class Parameters:
    """Class to store parameters like num_vg_idvd, Cox, etc"""
   
    def __init__(self, col_index, num_var_idvg, num_vd_idvg, W, Cox,\
                 num_var_idvd = None, num_vg_idvd = None, isbipolar_idvd = None,\
                     isbipolar_idvg = 1, Vds_param = None, decimals_Vg = 0, interp = 0,\
                     column = None, row = None):
        self.col_index = col_index
        self.num_var_idvg = num_var_idvg
        self.num_vd_idvg = num_vd_idvg
        self.num_var_idvd = num_var_idvd
        self.num_vg_idvd = num_vg_idvd
        self.W = W
        self.Cox = Cox
        self.isbipolar_idvd = isbipolar_idvd
        self.isbipolar_idvg = isbipolar_idvg
        self.Vds_param = Vds_param
        self.decimals_Vg = decimals_Vg
        self.interp = interp
        self.column = column
        self.row = row

class IdVg:
    """Class for Id-Vg data"""

    def __init__(self, inp, params, step_number):  # Argument is of type InputData
        """
        1) inp - input data object
        2) step_number - tells us which Vds step is to be considered
        """
        col_index = params.col_index
        num_var_idvg = params.num_var_idvg
        isbipolar_idvg = params.isbipolar_idvg
        W = params.W
        
        self.L = inp.channel_L  # Channel length (um)
        self.gateLeak = 0 # Flag indicating gate leakage
        self.compliance = 0 # Whether the device hit compliance for any of the currents
        
        if self.L < 1:
            self.L_label = str(int(self.L*1000))+r' nm'
        else:
            self.L_label = str(self.L)+r' $\mu$m'
        
        self.channel_length = inp.channel_length  # Channel length string
        self.num_vg = int(inp.id_vg.shape[0])  # Number of Vg points in Id-Vg
        
        if isbipolar_idvg:
            self.num_vg_forward = int(inp.id_vg.shape[0] / 2)  # Number of Vg points in                
        else:
            self.num_vg_forward = int(inp.id_vg.shape[0])  # Number of Vg points in      
            
        self.Id = inp.id_vg[0:self.num_vg, col_index[0]+step_number*num_var_idvg] / W  # Reading drain current (A/um)        
        self.Vg = inp.id_vg[0:self.num_vg, col_index[1]+step_number*num_var_idvg]  # Reading gate voltage (V)
        self.Ig = inp.id_vg[0:self.num_vg, col_index[6]+step_number*num_var_idvg]  # Reading gate current (A/um)        
        
        if np.max([np.max(np.abs(self.Id)),np.max(np.abs(self.Ig)),np.max(np.abs(self.Vg))]) > compliance_threshold:
            self.compliance = 1
            
        self.Id_forward = inp.id_vg[0:self.num_vg_forward, col_index[0]+step_number*num_var_idvg] / W #(A/um)
        # Parsing only forward sweep drain current (A/um)
        self.Vg_forward = inp.id_vg[0:self.num_vg_forward, col_index[1]+step_number*num_var_idvg] #(V)
        
        if params.interp:
            Vg_min = np.round(np.min(self.Vg),0)
            Vg_max = np.round(np.max(self.Vg),0)
            Vg_interp = np.arange(Vg_min, Vg_max+1) # Interpolate with dV step = 1
            Id_interp = np.interp(Vg_interp, self.Vg_forward, self.Id_forward)
            self.Id_forward = Id_interp #(A/um)
            self.Vg_forward = Vg_interp #(V)
            self.num_vg_forward = int(self.Id_forward.shape[0])  # Number of Vg points in      
        
        self.Id_max_idvg = np.amax(np.absolute(self.Id))*1e6 #uA/um
        # Parsing only forward sweep gate voltage (V)
        
        if params.Vds_param == None:
            self.Vds = np.round(inp.id_vg[0, col_index[2]+step_number*num_var_idvg], 1)  # Reading Vds from Id-Vg data (V)                
        else:
            self.Vds = params.Vds_param
        
        self.I_off = np.maximum(self.Id_forward[0], I_noise/W)  # Off current (A/um)        
        # (comparing with noise floor)
        self.on_off_ratio = []  # On-off ratio (no units)
        self.SS = []  # Subthreshold swing (mV/dec)
        self.mu_FE = []  # Field-effect mobility (cm2/V/s)
        self.gm = []  # Transconductance (S/um)
        self.dgm = []  # Transconductance derivative (S/V/um)
        self.Vt_const = []  # Constant current threshold voltage (V)
        self.Vt_lin = []  # Linear extrapolation threshold voltage (V)
        self.Vt_dgm = []  # Threshold voltage from max derivative of gm (V)
        self.Vov = []  # Overdrive voltage Vgs-Vt (V)
        self.goodIdVg = 1 # Flag to show whether Id-Vg curve is well-behaved or "good"


    def idvg_calc(self, params):
        """Calculation of parameters from Id-Vg"""
        
        Cox = params.Cox
        W = params.W
        self.on_off_ratio = self.Id_forward[-1] / self.I_off  # On-off ratio (no units)        
        # Off current is minimum(Id, noise floor=1e-12A)
        self.SS = np.diff(self.Vg_forward) / np.diff(np.log10(np.absolute(self.Id_forward.astype(float)) + np.spacing(1)) + np.spacing(1)) * 1e3  # Subthreshold slope (mV/dec)
        # Calculating field-effect mobility
        self.mu_FE = np.zeros((self.Vg_forward.shape[0] - 1))
        self.gm = np.diff(smooth(self.Id_forward, 5)) / np.diff(self.Vg_forward)  # Transconductance (S/um)
        self.dgm = np.diff(smooth(self.gm, 5)) / np.diff(self.Vg_forward[:-1])  # Transconductance derivative (S/V/um)
        self.dgm = smooth(self.dgm, 5)  # Smoothed Transconductance derivative (S/V/um)
        self.gm = smooth(self.gm)  # Smoothed Transconductance (S/um)
        self.mu_FE = self.L * self.gm / (Cox * self.Vds * 1e-9)  # Possible under-estimation? (cm2/V/s)
        
        # Calculating Constant Current Vt ##
        if np.amin(np.absolute(self.Id_forward)) < Id_thresh * W / self.L:
            self.Vt_const = np.amax(
                self.Vg_forward[self.Id_forward < Id_thresh * W / self.L])  # Constant current Vt (V)
        
        gmmax_index = np.where(self.gm==np.amax(self.gm))[0] #Index of maximum gm                                
        dgmmax_index = np.where(self.dgm[:-9] == np.amax(self.dgm[:-9]))  # Index of maximum dgm
        dgmmax_index, _ = find_peaks(self.dgm, height=0)  # Finds the first local max in dgm        
        
        if self.Id_max_idvg < 1e-1 and self.on_off_ratio < 1e2: #Setting gate leakage threshold at 1e-7
            self.goodIdVg = 0

        if (np.amax(np.absolute(self.Ig)) > 1e-7): #Setting gate leakage threshold at 1e-7
            self.gateLeak = 1
        
        if (dgmmax_index.shape[0] > 0):
            # Calculating Max Transconductance Derivative Vt ##        
            dgmmax_index = dgmmax_index[0]            
            self.Vt_dgm = self.Vg_forward[dgmmax_index]  # Threshold voltage from max derivative of gm (V)            
            Vt = self.Vt_dgm        
        if (gmmax_index.shape[0] > 0):
        # elif (gmmax_index.shape[0] > 0):
            # Calculating Linear Interpolation Vt ##
            gmmax_index = gmmax_index[0]        
            Vg_intercept = self.Vg_forward[gmmax_index] - self.Id_forward[gmmax_index] / np.amax(self.gm)
            Vt_lin_temp = Vg_intercept - self.Vds / 2  # Linear extrapolated Vt (V)        
            self.Vt_lin = find_nearest(self.Vg_forward, Vt_lin_temp) # Rounding the Vt_lin to the nearest element in Vg array        
            Vt = self.Vt_lin                
        elif (dgmmax_index.shape[0] > 0):
            # Calculating Max Transconductance Derivative Vt ##        
            dgmmax_index = dgmmax_index[0]            
            self.Vt_dgm = self.Vg_forward[dgmmax_index]  # Threshold voltage from max derivative of gm (V)            
            Vt = self.Vt_dgm        
        else:
            Vt = 0              
        self.Vov = self.Vg_forward - Vt #Overdrive voltage wrt linear Vt (V)        
        

    def plot_IdVg(self, ax, fig, lg='lin', anno = 'no'):
        """Plots Id-Vg characteristics. 
        1) lg == 'lin' for linear scale 
        2) lg == 'log' for log scale
        
        anno - whether to annotate or not"""

        if anno == 'yes':
            ax.annotate(r'$V_{DS}$ = ' + str(self.Vds) + ' V', xy=(.45, annotate_y1), xycoords='figure fraction')        
        ax.set(xlabel=vg_label, ylabel=id_label)
        if lg == 'lin':
            ax.locator_params(nbins=6, axis='y')            
        if lg == 'log':
            ax.set_yscale('log')
        ax.plot(self.Vg, self.Id * 1e6)  # Plotting Id (uA/um) vs Vgs (V)


    def plot_Ig(self, ax, fig, lg='log'):
        """Plots Id-Vg characteristics. 
        1) lg == 'lin' for linear scale 
        2) lg == 'log' for log scale
        
        anno - whether to annotate or not"""
        
        ax.set(xlabel=vg_label, ylabel=ig_label)
        if lg == 'lin':
            ax.locator_params(nbins=6, axis='y')            
        if lg == 'log':
            ax.set_yscale('log')
        ax.plot(self.Vg, np.abs(self.Ig))  # Plotting Id (uA/um) vs Vgs (V)        


    def plot_gm(self, ax, fig, legend = 'no'):
        """ Plots gm vs Vgs"""        
        
        ax.set(xlabel=vg_label, ylabel=gm_label)
        Vgs_temp = self.Vg_forward[:-1] #Since gm is derivative            
        ax.plot(Vgs_temp, self.gm)
        legends = [self.L_label[-1]]        
        if legend == 'yes':
            ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])
        return legends


class IdVd:
    """Class for Id-Vd data"""

    def __init__(self, inp, params):  # Argument is of type InputData
        col_index = params.col_index
        num_var_idvd = params.num_var_idvd
        num_vg_idvd = params.num_vg_idvd
        isbipolar_idvd = params.isbipolar_idvd
        W = params.W
        if isbipolar_idvd:  # If bipolar sweep, just take forward sweep, else dual sweep
            num_vd = int(inp.id_vd.shape[0])
        else:
            num_vd = int(inp.id_vd.shape[0])
        self.Vd = np.round(inp.id_vd[0:num_vd, col_index[5]], 2)  # Drain voltage (V)
        self.Id = inp.id_vd[0:num_vd,
                  col_index[3]:num_vg_idvd * num_var_idvd:num_var_idvd] / W  # Drain current (A/um)
        self.Vg = inp.id_vd[0, col_index[4]:num_vg_idvd * num_var_idvd:num_var_idvd]
        self.channel_length = inp.channel_length  # Channel length string
        self.L = inp.channel_L  # Channel length (um)
        
        if self.L < 1:
            self.L_label = str(int(self.L*1000))+r' nm'
        else:
            self.L_label = str(self.L)+r' $\mu$m'
        self.Id_max = np.amax(np.abs(
            inp.id_vd[:, col_index[3] + num_var_idvd * (num_vg_idvd - 1)] / W * 1e6))  # Maximum drain current (uA/um)


    def plot_IdVd(self, ax, fig, params, Vgs=None):
        """Plots Id-Vd characteristics in linear scale"""
        
        isbipolar_idvd = params.isbipolar_idvd        
        
        if Vgs is None:  # When plot_IdVd is called directly, not through TLM
            ax.set(xlabel=vd_label, ylabel=id_label)
            ax.annotate(r'$L_{ch}$ = ' + self.L_label, xy=(annotate_x, annotate_y1), xycoords='axes fraction')
            ax.plot(self.Vd, self.Id * 1e6)  # Plotting Id (uA/um) vs Vds (V) for alL Vgs            
            ax.locator_params(nbins=6, axis='y')            
            legends = [str(i) + ' V' for i in self.Vg]            
            ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.05])
            #Set tick frequency
            # set_ticks(ax, 5)
            if isbipolar_idvd == 1: # in case of bipolar Id-Vd, indicate origin
                ax.axhline(0, color='grey', ls='dotted', lw=0.6)
                ax.axvline(0, color='grey', ls='dotted', lw=0.6)  
        else:            
            i = np.nonzero(np.round(self.Vg, 1) == float(Vgs))                           
            i = i[0].item()  # Finding the index based on Vg
            ax.plot(self.Vd, self.Id[:, i] * 1e6)
            # Plotting Id (uA/um) vs Vds (V) for given Vgs                                                 
        
        
class TLM:
    """Class for TLM data"""

    def __init__(self, inp_data, inp_count, params):        
        self.params = params
        num_vd_idvg = self.params.num_vd_idvg
        self.column = self.params.column #Column of the TLM (1-9)
        self.row = self.params.row #Row of the TLM (1-15)
        self.data = inp_data  # Input data for all channel lengths
        self.count = inp_count  # Number of channel lengths        
        self.idvg = np.empty(shape=(self.count,num_vd_idvg), dtype=object)        
        self.idvd = np.empty(shape=(self.count,), dtype=object)  # Id-Vd object        
        self.L = np.empty(shape=(self.count,), dtype=float)  # Channel lengths in the TLM
        self.L_label = np.empty(shape=(self.count,), dtype=object)  # Label strings to use for plotting
        self.gateLeak = 0 #Flag to indicate gate leakage
        self.compliance = 0 #Flag to indicate whether device hit compliance        
        self.position = str(self.row)+'-'+str(self.column) #TLM position string
        
        goodIdVg_count = 0
        for i in np.arange(self.count):
            goodIdVg_flag = 1 # Flag to represent whether IdVg is well-behaved (1) or not (0)
            gateLeak_temp = self.gateLeak
            compliance_temp = self.compliance
            for j in np.arange(num_vd_idvg):
                self.idvg[goodIdVg_count,j] = IdVg(self.data[i],self.params, step_number = j)  # Id-Vg object
                self.idvg[goodIdVg_count,j].idvg_calc(self.params)  # Perform calculations on Id-Vg
                if self.idvg[goodIdVg_count,j].gateLeak:
                    self.gateLeak = 1
                if self.idvg[goodIdVg_count,j].compliance:
                    self.compliance = 1
                if self.idvg[goodIdVg_count,j].goodIdVg == 0:
                    goodIdVg_flag = 0
            if goodIdVg_flag == 0:
                self.gateLeak = gateLeak_temp
                self.compliance = compliance_temp
            elif goodIdVg_flag == 1:
                if self.data[goodIdVg_count].id_vd is None:
                    self.idvd[goodIdVg_count] = None
                    self.Vg_idvd = None  # Vg used in Id-Vd (V)
                else:
                    self.idvd[goodIdVg_count] = IdVd(self.data[i], self.params)  # Id-Vd object
                    self.Vg_idvd = self.idvd[0].Vg  # Vg used in Id-Vd (V)                
                self.L[goodIdVg_count] = self.idvg[goodIdVg_count,0].L
                if self.L[goodIdVg_count] < 1:
                    self.L_label[goodIdVg_count] = str(int(self.L[goodIdVg_count]*1000))+r' nm'
                else:
                    self.L_label[goodIdVg_count] = str(self.L[goodIdVg_count])+r' $\mu$m'
                goodIdVg_count = goodIdVg_count + 1
        
        self.count = goodIdVg_count                                            
        self.idvg = self.idvg[:self.count,:]
        self.idvd = self.idvd[:self.count]
        self.L = self.L[:self.count]
        self.L_label = self.L_label[:self.count]
        
        num_vg_forward = int(self.idvg[0,0].num_vg_forward)  # Number of Vg points in forward sweep of Id-Vg data
        # num_vg_forward = int(self.idvg[0].num_vg_forward)  # Number of Vg points in forward sweep of Id-Vg data
        num_vg = int(self.idvg[0,0].num_vg)  # Number of Vg points in Id-Vg data        
        # num_vg = int(self.idvg[0].num_vg)  # Number of Vg points in Id-Vg data        
        self.Vds_idvg = [self.idvg[0,0].Vds,  self.idvg[0,-1].Vds] # Vds forId-Vg sweep        
        self.Vov = np.zeros([num_vg_forward, self.count])  # Overdrive voltage (V)
        self.Vg_idvg = self.idvg[0,0].Vg  # Vg used in Id-Vg (V)                
        self.Vg_idvg_forward = self.idvg[0,0].Vg_forward  # Vg used in Id-Vg (V)                
        # Id from Id-Vg collected for all channel lengths(A/um)
        
        self.Id_adjusted = []  # Id adjusted after Vt alignment (A/um)
        self.mu_FE_adjusted = []  # mu_FE adjusted after Vt alignment (A/um)
        self.Rtot = []  # Total resistance (Ohm.um)
        self.Rc = []  # Contact resistance (Ohm.um)
        self.del_Rc = []  # Error in Contact resistance (Ohm.um)
        self.Rsh = []  # Sheet resistance (Ohm/sq)
        self.del_Rsh = []  # Error in Sheet resistance (Ohm/sq)
        self.mu_TLM = []  # Mobility extracted from TLM (cm2/V/s)
        self.n2D = []  # 2D carrier concentration (cm-2)
        self.TLM_fit = []  # Pandas dataframe to store regression data @ highest Vov
        self.R_squared = 0  # R-squared parameter for TLM fit
        self.goodTLM = 0 # Flag to indicate a well-behaved TLM
        self.goodRc = 0 #Flag to indicate good Rc extraction
        self.goodMobility = 0 #Flag to indicate good mu extraction 
                
        self.Id_idvg = np.zeros([num_vg, self.count,num_vd_idvg])        
        for i in np.arange(self.count):                         
            self.Id_idvg[:, i,0] = self.idvg[i,0].Id
            self.Id_idvg[:, i,-1] = self.idvg[i,-1].Id


    def tlm_calc(self, vd_step = 'low'):
        """Calculates Rtot,Rc,Rsh,mu_TLM
        vd_step = 'low' for low_vds
        vd_step = 'hi' for high_vds"""        
        
        Cox = self.params.Cox
        if vd_step == 'low':
            j = 0
        elif vd_step == 'hi':
            j = -1
        # Finding Vgs-Vt overlap and re-arranging
        
        for i in np.arange(self.count):
            self.Vov[:, i] = self.idvg[i,j].Vov              
                              
        self.Vov = np.round(self.Vov, self.params.decimals_Vg)        
        Vov_hi = np.amin(self.Vov.max(0))  # Higher and lower limits of Vov (or Vgs-Vt)
        Vov_lo = np.amax(self.Vov.min(0))
        index = np.logical_and(self.Vov < Vov_hi, self.Vov > Vov_lo)
        
        dV = np.mean(np.diff(self.Vov[:,0]))        
        Vov_interp = np.arange(Vov_lo, Vov_hi, dV)
       
        # Window of Vgs-Vt overlap across multiple Lchannel's
        self.Vov = Vov_interp
        self.Id_adjusted = np.zeros([self.Vov.shape[0], self.count]) #(A/um)
        self.mu_FE_adjusted = np.zeros([self.Vov.shape[0]-1, self.count]) #(cm2/V/s)
        self.gm_adjusted = np.zeros([self.Vov.shape[0]-1, self.count]) #(S/um)
        
        if self.Vov.shape[0] > 0 and self.count > 1 and self.gateLeak == 0 and self.compliance == 0:            
            self.goodTLM = 1
        
        # Adjusting Id to align the Vgs-Vt for various channels                        
        
        if self.goodTLM:
            for i in np.arange(self.count):                            
                self.Id_adjusted[:, i] = np.interp(self.Vov, self.idvg[i,j].Vov, self.idvg[i,j].Id_forward) #(A/um)       
                self.mu_FE_adjusted[:, i] = np.interp(self.Vov[:-1], self.idvg[i,j].Vov[:-1], self.idvg[i,j].mu_FE) #(cm2/V/s)                                
                self.gm_adjusted[:, i] = np.interp(self.Vov[:-1], self.idvg[i,j].Vov[:-1], self.idvg[i,j].gm) #(S/um)                

            # Calculating Rc, del_Rc, Rsh, del_Rsh
            self.Rtot = self.Vds_idvg[j] / self.Id_adjusted # Total resistance (Ohm.um)            
            self.Rc = np.zeros(self.Vov.shape[0])           #(Ohm.um)
            self.del_Rc = np.zeros(self.Vov.shape[0])       #(Ohm.um)
            self.Rsh = np.zeros(self.Vov.shape[0])          #(Ohm.um)
            self.del_Rsh = np.zeros(self.Vov.shape[0])      #(Ohm.um)
        
            # TLM fitting
                    
            for i in np.arange(self.Vov.shape[0]):
                X = sm.add_constant(self.L.conj().transpose())  # Adding constant term to the model                
                model = sm.OLS(self.Rtot[i, :].conj().transpose(), X)                
                results = model.fit()
                self.Rc[i] = results.params[0] / 2              # Contact resistance (Ohm.um)
                self.del_Rc[i] = results.bse[0] / 2             # Error in Contact resistance (Ohm.um)
                if ~(np.isfinite(self.del_Rc[i])):
                    self.del_Rc[i] = 0                              
                self.Rsh[i] = results.params[1]                 # Sheet resistance (Ohm/sq)                
                self.del_Rsh[i] = results.bse[1]                # Error in Sheet resistance (Ohm/sq)
                if ~(np.isfinite(self.del_Rsh[i])):
                    self.del_Rsh[i] = 0              
            self.R_squared = results.rsquared
            d = {rtot_label: self.Rtot[i, :].conj().transpose() / 1e3,
                 lch_label: self.L.conj().transpose()}  # Creating dictionary of Rtot vs L
            self.TLM_fit = pd.DataFrame(data=d)  # Converting dictionary to Pandas dataframe        
            self.n2D = Cox * 1e-9 * (self.Vov + np.spacing(1)) / q  # 2D carrier concentration (cm-2)        
            # np.spacing(1) is used to avoid divide-by-zero error in next step
            self.mu_TLM = 1 / ((self.Rsh * self.n2D) * q + np.spacing(1))  # Mobility extracted from TLM (cm^2/V/s)
            self.del_mu_TLM = self.del_Rsh / (self.Rsh + np.spacing(1)) * self.mu_TLM  # Error in mobility extracted from TLM (cm^2/V/s)
            if (self.Rc[-1] > 0) and (self.Rc[-1] > 2*np.abs(self.del_Rc[-1])):
                self.goodRc = 1
            if (self.mu_TLM[-1] > 0) and (self.mu_TLM[-1] > 2*np.abs(self.del_mu_TLM[-1])):
                self.goodMobility = 1
            if np.max(self.n2D) <= n2D_LL:
                self.goodTLM = 0


    def plot_TLM_fit(self, ax, fig):
        """Plots TLM fit - regression plot"""            
        
        if self.goodTLM == 0:
            return
                        
        plt.gca().set_prop_cycle(None)
        ax.set_xlim(left=0, right=round(self.L[-1]*1.1/0.2)*0.2)  # Sets lower x-limit to 0
        sns.regplot(x=lch_label, y=rtot_label, data=self.TLM_fit,\
                        ci=95, truncate=False, ax=ax)         
        # truncate = False plots the regression line to fill x-limits
        start, end = ax.get_ylim()
        ax.set_ylim(bottom=np.minimum(0, start))  # Sets lower y-limit to <=0
        ax.set_ylim(top=np.maximum(3*start, end))  # Sets upper y-limit        
        Rc_anno = np.round(self.Rc[-1] / 1e3, 1)
        delRc_anno = np.round(self.del_Rc[-1] / 1e3, 1)
        ax.annotate(r'$R_{C}$ = ' + str(Rc_anno) + r' $\pm$ ' + str(delRc_anno) + r' k$\Omega$.$\mu$m',
                    xy=(annotate_x, annotate_y1), xycoords='axes fraction')
        n2D_anno = self.n2D[-1] / 1e12
        if n2D_anno < 10:
            ax.annotate(r'at $n_{2D}$ = ' + str(np.round(n2D_anno, 1)) + r' $\times$ $10^{12} cm^{-2}$',
                    xy=(annotate_x, annotate_y22), xycoords='axes fraction', fontsize=8.5)      
        else:
            ax.annotate(r'at $n_{2D}$ = ' + str(np.round(n2D_anno/10, 1)) + r' $\times$ $10^{13} cm^{-2}$',
                    xy=(annotate_x, annotate_y22), xycoords='axes fraction', fontsize=8.5)      
        ax.annotate(r'R-squared = ' + str(np.round(self.R_squared, 3)),
                    xy=(.6, .1), xycoords='axes fraction', fontsize=8.5)
        # set_ticks(ax,5)
        # set_ticks(ax,5,'x',0.2)


    def plot_mu_FE(self, ax, fig, channel = 0, flag=None, legend = 'no'):
        """Plots mu_FE vs Vov or n2D for the specified channel
            1) flag == None - plot vs Vov
            2) flag == 'n2D'- plot vs n2D
            
            channel = 0: plot all channels
            channel = -1: plot longest channel
            channel = n: n-th channel"""        
        # Plotting mu_FE (cm2/V/s) vs Vov(V)                
        if flag is None:
            ax.set(xlabel=vov_label, ylabel=mu_fe_label)
            Vov_temp = self.Vov[:-1] #Since mu_FE is derivative
            if channel == 0:
                """plot mu_FE for all channels"""         
                for i in np.arange(self.count):    
                    ax.plot(Vov_temp[Vov_temp >= 0], self.mu_FE_adjusted[Vov_temp >= 0, i])
                legends = [i for i in self.L_label]
            elif channel == -1:
                """plot mu_FE for longest channels"""         
                ax.plot(Vov_temp[Vov_temp >= 0], self.mu_FE_adjusted[Vov_temp >= 0, -1])
                legends = [self.L_label[-1]]
            elif channel > 0:
                """plot mu_FE for the specified channel"""
                ax.plot(Vov_temp[Vov_temp >= 0], self.mu_FE_adjusted[Vov_temp >= 0, channel-1])
                legends = [self.L_label[channel-1]]
        # Plotting mu_FE (cm2/V/s) vs n2D(V)            
        elif flag == 'n2D':
            n2D_temp = self.n2D[:-1] #Since mu_FE is derivative
            ax.set(xlabel=n2d_label, ylabel=mu_fe_label)
            if channel == 0:
                """plot mu_FE for all channels"""
                for i in np.arange(self.count):
                    ax.plot(n2D_temp[n2D_temp >= n2D_LL] / 1e12, self.mu_FE_adjusted[n2D_temp >= n2D_LL, i])
                legends = [i for i in self.L_label]
            elif channel == -1:
                """plot mu_FE for longest channels"""         
                ax.plot(n2D_temp[n2D_temp >= n2D_LL] / 1e12, self.mu_FE_adjusted[n2D_temp >= n2D_LL, -1])
                legends = [self.L_label[-1]]
            elif channel > 0:
                """plot mu_FE for the specified channel"""
                ax.plot(n2D_temp[n2D_temp >= n2D_LL] / 1e12, self.mu_FE_adjusted[n2D_temp >= n2D_LL, channel-1])
                legends = [self.L_label[channel-1]]
        if legend == 'yes':
            ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])
        return legends


    def plot_mu_TLM(self, ax, fig, flag=None):
        """Plots mu_TLM vs Vov or n2D
            1) flag == None - plot vs Vov
            2) flag == 'n2D'- plot vs n2D"""        
    
        if self.goodTLM == 0:
            return 
        
        if flag is None:
            ax.set(xlabel=vov_label, ylabel=rc_label)
            ax.errorbar(self.Vov[self.Vov >= 0], self.mu_TLM[self.Vov >= 0],
                        self.del_mu_TLM[self.Vov >= 0])
            # Plotting mu_TLM (cm2/V/s) vs Vov(V)
        elif flag == 'n2D':
            ax.set(xlabel=n2d_label, ylabel=mu_eff_label)
            ax.errorbar(self.n2D[self.n2D >= n2D_LL] / 1e12, self.mu_TLM[self.n2D >= n2D_LL],
                        self.del_mu_TLM[self.n2D >= n2D_LL])                
    
        mu_TLM_cut = self.mu_TLM[self.n2D >= n2D_LL]
        del_mu_TLM_cut = self.del_mu_TLM[self.n2D >= n2D_LL]
        top_limit = np.maximum(2*mu_TLM_cut[0]\
                                            ,mu_TLM_cut[0]+del_mu_TLM_cut[0])
        ax.set_ylim(bottom=0,top=top_limit)  # Sets lower y-limit to 0
        x_LL = np.maximum(n2D_LL,self.n2D[0])
        ax.set_xlim(left=x_LL/1e12)
        mu_TLM_anno = int(np.round(self.mu_TLM[-1]))
        del_mu_TLM_anno = int(np.round(self.del_mu_TLM[-1]))
        ax.annotate(r'$\mu_{eff}$ = ' + str(mu_TLM_anno) + r' $\pm$ ' +\
                    str(del_mu_TLM_anno) + r' $cm^2$$V^{-1}$$s^{-1}$',
                    xy=(annotate_x, annotate_y1), xycoords='axes fraction')
        n2D_anno = self.n2D[-1] / 1e12
        if n2D_anno < 10:
            ax.annotate(r'at $n_{2D}$ = ' + str(np.round(n2D_anno, 1)) + r' $\times$ $10^{12} cm^{-2}$',
                    xy=(annotate_x, annotate_y22), xycoords='axes fraction', fontsize=8.5)      
        else:
            ax.annotate(r'at $n_{2D}$ = ' + str(np.round(n2D_anno/10, 1)) + r' $\times$ $10^{13} cm^{-2}$',
                    xy=(annotate_x, annotate_y22), xycoords='axes fraction', fontsize=8.5)                  
        # set_ticks(ax,5)
        

    def plot_mu_compare(self, ax, fig, channel = 0, flag=None, legend = 'no'):
        """Compare mu_TLM and mu_FE"""
        legends = self.plot_mu_FE(ax, fig, channel, flag, legend = 'no')
        self.plot_mu_TLM(ax, fig, flag)
        legends = legends + [r'$\mu_{eff}$']
        if legend == 'yes':            
            ax.legend(legends, loc='upper right',bbox_to_anchor=[0.95,0.95])


    def plot_gm(self, ax, fig, channel = 0, legend = 'no'):
        """ Plot gm vs Vgs
            channel = 0: plot all channels
            channel = -1: plot longest channel
            channel = n: n-th channel"""    
            
        if channel == 0:
            """plot gm for all channels"""         
            for i in np.arange(self.count):    
                self.idvg[i, 0].plot_gm(ax, fig)
            legends = [i for i in self.L_label]
        elif channel == -1:
            """plot gm for longest channels"""         
            self.idvg[-1, 0].plot_gm(ax, fig)
            legends = [self.L_label[-1]]
        elif channel > 0:
            """plot gm for the specified channel"""
            self.idvg[channel-1, 0].plot_gm(ax, fig)
            legends = [self.L_label[channel-1]]
        if legend == 'yes':
            ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])
        return legends


    def plot_gm_aligned(self, ax, fig, channel = 0, flag=None, legend = 'no'):
        """ Plot gm vs n2D/Vov, aligned wrt Vth
            1) flag == None - plot vs Vov
            2) flag == 'n2D'- plot vs n2D
            
            channel = 0: plot all channels
            channel = -1: plot longest channel
            channel = n: n-th channel"""        
        # Plotting gm (cm2/V/s) vs Vov(V)                
        if flag is None:
            ax.set(xlabel=vov_label, ylabel=gm_label)
            Vov_temp = self.Vov[:-1] #Since gm is derivative
            if channel == 0:
                """plot mu_FE for all channels"""         
                for i in np.arange(self.count):    
                    ax.plot(Vov_temp[Vov_temp >= 0], self.gm_adjusted[Vov_temp >= 0, i])
                legends = [i for i in self.L_label]
            elif channel == -1:
                """plot gm for longest channels"""         
                ax.plot(Vov_temp[Vov_temp >= 0], self.gm_adjusted[Vov_temp >= 0, -1])
                legends = [self.L_label[-1]]
            elif channel > 0:
                """plot gm for the specified channel"""
                ax.plot(Vov_temp[Vov_temp >= 0], self.gm_adjusted[Vov_temp >= 0, channel-1])
                legends = [self.L_label[channel-1]]
        # Plotting gm (cm2/V/s) vs n2D(V)            
        elif flag == 'n2D':
            n2D_temp = self.n2D[:-1] #Since gm is derivative
            ax.set(xlabel=n2d_label, ylabel=gm_label)
            if channel == 0:
                """plot gm for all channels"""
                for i in np.arange(self.count):
                    ax.plot(n2D_temp[n2D_temp >= n2D_LL] / 1e12, self.gm_adjusted[n2D_temp >= n2D_LL, i])
                legends = [i for i in self.L_label]
            elif channel == -1:
                """plot gm for longest channels"""         
                ax.plot(n2D_temp[n2D_temp >= n2D_LL] / 1e12, self.gm_adjusted[n2D_temp >= n2D_LL, -1])
                legends = [self.L_label[-1]]
            elif channel > 0:
                """plot gm for the specified channel"""
                ax.plot(n2D_temp[n2D_temp >= n2D_LL] / 1e12, self.gm_adjusted[n2D_temp >= n2D_LL, channel-1])
                legends = [self.L_label[channel-1]]
        if legend == 'yes':
            ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])
        return legends


    def get_Rc(self, n2D = 1e13):
        """Returns Rc at specified n2D"""        
        
        Rc_interp = np.interp(n2D, self.n2D, self.Rc)
        del_Rc_interp = np.interp(n2D, self.n2D, self.del_Rc)
        goodRc = 0
        if Rc_interp > 0 and Rc_interp < 1e4 and np.abs(Rc_interp) > np.abs(del_Rc_interp):
            goodRc = 1
        return Rc_interp, del_Rc_interp, self.R_squared, goodRc, rc_cdf_label
    
    
    def get_mu_TLM(self, n2D = 1e13):
        """Returns mu_TLM at specified n2D"""        
        
        mu_TLM_interp = np.interp(n2D, self.n2D, self.mu_TLM)
        del_mu_TLM_interp = np.interp(n2D, self.n2D, self.del_mu_TLM)
        goodmu_TLM = 0
        if mu_TLM_interp > 0 and np.abs(mu_TLM_interp) > np.abs(del_mu_TLM_interp):
            goodmu_TLM = 1
        return mu_TLM_interp, del_mu_TLM_interp, self.R_squared, goodmu_TLM, mu_eff_label


    def plot_Rtot(self, ax, fig, n2D = 1e13, anno = 'yes'):
        """Plots Rtot at specified n2D"""        
        
        if anno == 'yes':
            n2D_anno = n2D / 1e12
            if n2D_anno < 10:
                ax.annotate(r'at $n_{2D}$ = ' + str(np.round(n2D_anno, 1)) + r' $\times$ $10^{12} cm^{-2}$',
                    xy=(annotate_x, annotate_y22), xycoords='axes fraction', fontsize=8.5)      
            else:
                ax.annotate(r'at $n_{2D}$ = ' + str(np.round(n2D_anno/10, 1)) + r' $\times$ $10^{13} cm^{-2}$',
                    xy=(annotate_x, annotate_y22), xycoords='axes fraction', fontsize=8.5)                          
        for j in np.arange(self.count):            
            Rtot_interp = np.interp(n2D, self.n2D, self.Rtot[:,j])
            ax.plot(self.L[j], np.abs(Rtot_interp/1e3),'o',color='black') #Plotting Rtot(kohm.um) vs Lch (um)
        pass
    
    
    def plot_Rc(self, ax, fig, flag=None):
        """Plots Rc vs Vov or n2D
            1) flag == None - plot vs Vov
            2) flag == 'n2D'- plot vs n2D"""    
        
        if self.goodTLM == 0:
            return
        
        if flag is None:
            ax.set(xlabel=vov_label, ylabel=rc_label)
            ax.errorbar(self.Vov[self.Vov >= 0], self.Rc[self.Vov >= 0] / 1e3,
                        self.del_Rc[self.Vov >= 0] / 1e3)
            # Plotting Rc (kOhm.um) vs Vov(V)
        elif flag == 'n2D':
            ax.set(xlabel=n2d_label, ylabel=rc_label)
            ax.errorbar(self.n2D[self.n2D >= n2D_LL] / 1e12, self.Rc[self.n2D >= n2D_LL] / 1e3,
                        self.del_Rc[self.n2D >= n2D_LL] / 1e3)
            # Plotting Rc (kOhm.um) vs n2D(cm-2)          
            # Setting lower y-limit to min(0,Rc-del_Rc)                          
            try:
                ax.set_ylim(bottom=(np.minimum(np.amin(self.Rc[self.n2D >= n2D_LL]\
                                 - self.del_Rc[self.n2D >= n2D_LL]) / 1e3, 0)))            
            except ValueError:
                pass
            # Setting lower x-limit to max(n2D_LL,n2D[0])              
            ax.set_xlim(left=np.maximum(n2D_LL,self.n2D[0])/1e12)
        # set_ticks(ax,5)


    def plot_Rsh(self, ax, fig, flag=None):
        """Plots Rsh vs Vov or n2D
            1) flag == None - plot vs Vov
            2) flag == 'n2D'- plot vs n2D"""
        
        if self.goodTLM == 0:
            return
        
        if flag is None:
            ax.set(xlabel=vov_label, ylabel=rsh_label)
            ax.errorbar(self.Vov[self.Vov >= 0], self.Rsh[self.Vov >= 0] / 1e3,
                        self.del_Rsh[self.Vov >= 0] / 1e3)
            # Plotting Rc (kOhm.um) vs Vov(V)
        elif flag == 'n2D':
            ax.set(xlabel=n2d_label, ylabel=rsh_label)
            ax.errorbar(self.n2D[self.n2D >= n2D_LL] / 1e12, self.Rsh[self.n2D >= n2D_LL] / 1e3,
                        self.del_Rsh[self.n2D >= n2D_LL] / 1e3)
            ax.set_ylim(bottom=(np.minimum(np.amin(self.Rsh[self.n2D >= n2D_LL]\
                                 - self.del_Rsh[self.n2D >= n2D_LL]) / 1e3, 0)))            
            ax.set_xlim(left=np.maximum(n2D_LL,self.n2D[0])/1e12)
            Rsh_anno = int(np.round(self.Rsh[-1]/1000))
            del_Rsh_anno = int(np.round(self.del_Rsh[-1]/1000))
            ax.annotate(r'$R_{SH}$ = ' + str(Rsh_anno) + r'$\pm$' +\
                    str(del_Rsh_anno) + r' k$\Omega$.sq.',
                    xy=(.25, annotate_y1), xycoords='axes fraction')
            ax.annotate(r'at $n_{2D}$ = ' + str(np.round(self.n2D[-1] / 1e12, 1))\
                    + r' $\times$ $10^{12} cm^{-2}$',
                    xy=(.25, annotate_y22), xycoords='axes fraction', fontsize=8.5)


    def plot_IdVg(self, ax, fig, lg='lin', vd_step='low'):
        """Plots Id-Vg for all channel lengths in log scale
            1a)lg == None for linear scale 
            1b)lg == 'log' for log scale
            2a)vd_step == 'low' - low Vds
            2b)vd_step == 'hi' - hi Vds"""        
    
        ax.set(xlabel=vg_label, ylabel=id_label)
        if lg == 'log':
            ax.set_yscale('log')
            
        if vd_step == 'low':
            j = 0
            # ax.annotate(r'$V_{DS}$ = ' + str(self.Vds_idvg[j]) + ' V', xy=(.65, .65), xycoords='axes fraction')
            ax.annotate(r'$V_{DS}$ = ' + str(self.Vds_idvg[j]) + ' V', xy=(annotate_x, annotate_y1), xycoords='axes fraction')
            ax.plot(self.Vg_idvg, self.Id_idvg[:,:,j] * 1e6)  # Plotting Id (uA/um) vs Vgs (V)            
        elif vd_step == 'hi':            
            j = -1
            # ax.annotate(r'$V_{DS}$ = ' + str(self.Vds_idvg[j]) + ' V', xy=(.65, .65), xycoords='axes fraction')
            ax.annotate(r'$V_{DS}$ = ' + str(self.Vds_idvg[j]) + ' V', xy=(annotate_x, annotate_y1), xycoords='axes fraction')
            ax.plot(self.Vg_idvg, self.Id_idvg[:,:,j] * 1e6)  # Plotting Id (uA/um) vs Vgs (V)
        # ax.set_ylim(1e-1, 1e3)
        # ax.set_ylim(1e-6, 1e2)
        legends = [i for i in self.L_label]
        ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])


    def plot_IdVg_multVds(self, ax, fig, lg='lin',channel=0, anno='no'):
        """
        channel = 0-5 corresponds to CH1 - CH6
        anno = 'no' - no annotation; 'yes' - annotate
        """
        
        if anno == 'yes':
            ax.annotate(r'$L_{CH}$ = ' + self.idvg[channel,0].L_label, xy=(annotate_x, annotate_y1), xycoords='axes fraction')        
        num_vd_idvg = self.params.num_vd_idvg        
        self.idvg[channel,0].plot_IdVg(ax, fig, lg)  # Plotting Id-Vg for all channel lengths in log scale
        if num_vd_idvg > 1:        
            self.idvg[channel,-1].plot_IdVg(ax, fig, lg)  # Plotting Id-Vg for all channel lengths in log scale
            legends = [r'$V_{DS}$ = ' + str(self.Vds_idvg[0]) + ' V', r'$V_{DS}$ = ' + str(self.Vds_idvg[-1]) + ' V']
        else:
            legends = [r'$V_{DS}$ = ' + str(self.Vds_idvg[0]) + ' V']
        # ax.set_ylim(1e-6, 1e11)
        ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])
        set_x_ticks(ax,5)
        

    def plot_IdVg_Vov(self, ax, fig, lg='lin'):
        """Plots Id-Vg for all channel lengths in log scale
            1)lg == None for linear scale 
            2)lg == 'log' for log scale"""
    
        if self.goodTLM == 0:
            return    
        
        ax.annotate(r'$V_{DS}$ = ' + str(self.Vds_idvg) + ' V', xy=(.15, annotate_y1), xycoords='axes fraction')
        ax.set(xlabel=vov_label, ylabel=id_label)
        
        if lg == 'log':
            ax.set_yscale('log')
        ax.plot(self.Vov, self.Id_adjusted * 1e6)  # Plotting Id (uA/um) vs Vgs (V)
        legends = [i for i in self.L_label]
        ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])


    def plot_IdVd(self, ax, fig, Vgs):
        """Plots Id-Vd of for all channel lengths at a given Vgs""" 
        
        isbipolar_idvd = self.params.isbipolar_idvd

        for i in np.arange(self.count):
            self.idvd[i].plot_IdVd(ax, fig, params = self.params, Vgs = Vgs)        
        legends = [i for i in self.L_label]
        idmax_anno = int(np.round(np.amax(np.abs(self.idvd[0].Id*1e6))))
        ax.annotate(r'$V_{GS}$ = ' + str(Vgs) + ' V', xy=(0.03, annotate_y1), xycoords='axes fraction')        
        ax.annotate(r'$I_{D,max}$ = ' + str(idmax_anno)\
         + r' $\mu A / \mu m$', xy=(0.03, annotate_y22), xycoords='axes fraction')
        ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])
        ax.set(xlabel=vd_label, ylabel=id_label)
        ax.locator_params(nbins=6, axis='y')        
        #Set tick frequency
        # set_ticks(ax,5)
            
        if isbipolar_idvd == 1:
            ax.axhline(0, color='grey', ls='dotted', lw=0.6)
            ax.axvline(0, color='grey', ls='dotted', lw=0.6)


    def plot_Ig(self, ax, fig, lg='log'):
        """Plot Ig for all channels in the TLM"""
        pass
    
    
    def summary_plot(self, ax, fig, device_name):
        """
        Summary plot for the TLM:
        1) Id-Vg for all channel lengths at low-Vds
        2) Id-Vg for all channel lengths at hi-Vds
        3) Id-Vg for lowest channel length at both Vds
        4) TLM fit plot
        5) Rc-Vg plot
        6) mu_TLM plot
        """
        self.plot_IdVg(ax[0], fig, lg = 'log', vd_step = 'low')  # Plotting Id-Vg for all channel lengths in log scale
        self.plot_IdVg(ax[1], fig, lg = 'log', vd_step = 'hi')  # Plotting Id-Vg for all channel lengths in log scale
        self.plot_IdVg_multVds(ax[2], fig, 'log', channel = 0)  # Plotting Id-Vg for 100nm in log 
        if self.goodTLM:
            self.plot_TLM_fit(ax[3], fig)  # Plotting TLM fit for the last Vov
            self.plot_Rc(ax[4], fig, flag = 'n2D')  # Plotting Rc vs n2D
            self.plot_mu_TLM(ax[5], fig, flag = 'n2D')  # Plotting mu_eff vs n2D
        # self.plot_mu_FE(ax[5], fig, flag = 'n2D')  # Plotting mu_eff vs n2D        
        fig.suptitle(device_name)
        # tlm_set[device_count].plot_IdVg(axs_reshape[device_count],fig,'log') #Plotting Id-Vg for all channel lengths in log scale
        # axs_reshape[device_count].set_title(device[i])
        # tlm_set[device_count].plot_Rc(axs2_reshape[device_count],fig2,'n2D') #Plotting Rc vs n2Ds
        # axs2_reshape[device_count].set_title(device[i])
        # tlm_set[device_count].plot_mu_TLM(axs3_reshape[device_count],fig3,'n2D') #Plotting mu_TLM vs n2D
        # axs3_reshape[device_count].set_title(device[i])
        # tlm_set[device_count].plot_TLM_fit(axs4_reshape[device_count],fig4) #Plotting TLM fits
        # axs4_reshape[device_count].set_title(device[i])
        # device_count = device_count + 1            