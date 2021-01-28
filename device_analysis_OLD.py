"""
Created on Tue Oct 22 18:00:05 2019

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

# Choose plotting colors here
pal = sns.color_palette("bright")
# sns.palplot(sns.color_palette("bright", 8))
custom_colors = pal.as_hex()
# custom_colors = custom_colors[:1]+custom_colors[2:5]+custom_colors[6:7]+custom_colors[-1:]
custom_colors = custom_colors[:5]+custom_colors[6:7]+custom_colors[-1:]
# custom_colors = ["DarkBlue", "Crimson", "DarkGreen", "DarkOrchid", "Orange", "DeepPink"]
# You can use any HTML colors in the above line. For HTML colors, see this: https://www.w3schools.com/colors/colors_names.asp
sns.set_palette(custom_colors)

# Setting current directory and target directory
current_dir = os.getcwd()
# user_folder = os.getenv('USERPROFILE');
# user_folder = user_folder + "\\Google Drive\\Research\\Projects"
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
#dir_path = "\\Pd Interlayer contacts\\Janis\\2020-02-27-IA3-ALOX\\"; Vg_idvd_plt = 30
#dir_path = "\\InGaAs contacts\\Cascade\\2020-03-03-MF1"; Vg_idvd_plt = 70
dir_path = "\\Pd Interlayer contacts\\Janis\\2020-03-02-IB2\\"; Vg_idvd_plt = 20
# dir_path = "\\Pd Interlayer contacts\\Janis\\2020-08-30-IB13\\"; Vg_idvd_plt = 30
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

# Reading device parameters - file paths, column numbers, polarity of sweep, etc.
device_params = Path('device_parameters.xlsx')  # Device parameters file - what sweeps, column indices, etc.
device_params = xlrd.open_workbook(device_params, logfile=open(os.devnull, 'w'))
device_params = np.asarray(pd.read_excel(device_params, engine='xlrd'))  # Reading Id-Vg data
device = device_params[:, 2]  # List of devices
file_prefix = device_params[:, 3]  # List of file prefixes
isfolder_organized = device_params[:, 4]  # Whether the device data is organized into folders
col_index = device_params[np.ix_([0], [5, 6, 7, 8, 9, 10])]  # Column indices for various values
col_index = col_index[0]  # Converting 2-D array to 1-D
#1,2,3 - Ids, Vgs, Vds for Id-Vg
#4,5,6 - Ids, Vgs, Vds for Id-Vd
num_vd_idvg = device_params[0, 16]  # Number of Vd steps in Id-Vg sweep
num_var_idvg = device_params[0, 17]  # Number of column variables in a single Id-Vg sweep
num_vg_idvd = device_params[0, 11]  # Number of Vg steps in Id-Vd sweep
isbipolar_idvd = device_params[0, 12]  # Whether Id-Vd sweep is bipolar
num_var_idvd = device_params[0, 13]  # Number of column variables in a single Id-Vd sweep

W = device_params[0, 14]  # Channel width (um)
tox = device_params[0, 15]  # Thickness of the gate oxide (nm)

# Constants
q = 1.6e-19  # Charge of an electron (C)
Cox_30nm = 116  # 30nm SiO2 capacitance (Kirby ACS Nano) (nF/cm2)
Cox = Cox_30nm * 30 / tox  # 100nm SiO2 capacitance (nF/cm2)
Id_thresh = 2e-8  # For constant Vt calculation (A)
I_noise = 1e-12  # Noise floor for Id (A)
n2D_LL = 1e12  # n2D lower limit for plotting (cm-2)

# Plot labels
vd_label = r'Drain Voltage,  $V_{DS}$  [V]'
vg_label = r'Gate Voltage,  $V_{GS}$  [V]'
vov_label = r'Overdrive Voltage,  $V_{ov}$  [V]'
n2d_label = r'Carrier Concentration,  $n_{2D}$  [$10^{12}$ $cm^{-2}$]'
id_label = r'Drain Current,  $I_D$  [$\mu A$/$\mu m$]'
rc_label = r'Contact Resistance,  $R_{C}$  [k$\Omega.\mu m$]'
mu_eff_label = r'Effective Mobility,  $\mu_{eff}$  [$cm^2$$V^{-1}$$s^{-1}$]'
rsh_label = r'Sheet Resistance, $R_{sh}$ [k$\Omega$.sq.]'
rtot_label = r'Total Resistance, $R_{tot}$ [k$\Omega.\mu $m]'
lch_label = r'Channel Length,  $L_{ch}$  [$\mu$m]'
annotate_x = 0.05
annotate_y1 = 0.85
annotate_y2 = 0.75
annotate_y22 = 0.725

device_index = np.arange(device.shape[0])  # Device index

# Helper functions
def smooth(a, WSZ=7):
    """ a: NumPy 1-D array containing the data to be smoothed
    WSZ: smoothing window size needs, which must be odd number,
    as in the original MATLAB implementation"""
    out0 = np.convolve(a, np.ones(WSZ, dtype=int), 'valid') / WSZ
    r = np.arange(1, WSZ - 1, 2)
    start = np.cumsum(a[:WSZ - 1])[::2] / r
    stop = (np.cumsum(a[:-WSZ:-1])[::2] / r)[::-1]
    return np.concatenate((start, out0, stop))

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

class InputData:
    """Class for input data"""

    def __init__(self, idvg, idvd, index):
        self.id_vg = idvg  # Assigning read Id-Vg data
        self.channel_length = channel_length[index]  # Channel length (string)
        self.channel_L = channel_L[index]  # Channel length (value)
        self.id_vd = idvd  # Assigning read Id-Vd data


class IdVg:
    """Class for Id-Vg data"""

    def __init__(self, inp, step_number):  # Argument is of type InputData
        """
        1) inp - input data object
        2) step_number - tells us which Vds step is to be considered
        """
        self.L = inp.channel_L  # Channel length (um)
        if self.L < 1:
            self.L_label = str(int(self.L*1000))+r' nm'
        else:
            self.L_label = str(self.L)+r' $\mu$m'
        self.channel_length = inp.channel_length  # Channel length string
        self.num_vg = int(inp.id_vg.shape[0])  # Number of Vg points in Id-Vg
        self.num_vg_forward = int(inp.id_vg.shape[0] / 2)  # Number of Vg points in                
        self.Id = inp.id_vg[0:self.num_vg, col_index[0]+step_number*num_var_idvg] / W  # Reading drain current (A/um)
        self.Vg = inp.id_vg[0:self.num_vg, col_index[1]+step_number*num_var_idvg]  # Reading gate voltage (V)
        self.Id_forward = inp.id_vg[0:self.num_vg_forward, col_index[0]+step_number*num_var_idvg] / W
        # Parsing only forward sweep drain current (A/um)
        self.Vg_forward = inp.id_vg[0:self.num_vg_forward, col_index[1]+step_number*num_var_idvg]
        self.Id_max_idvg = np.amax(np.absolute(self.Id))*1e6 #uA/um
        # Parsing only forward sweep gate voltage (V)
        self.Vds = np.round(inp.id_vg[0, col_index[2]+step_number*num_var_idvg], 1)  # Reading Vds from Id-Vg data (V)                
        # self.Vds = 0.5
        self.I_off = np.maximum(self.Id[0], I_noise)  # Finding off current
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

    def idvg_calc(self):
        """Calculation of parameters from Id-Vg"""
        self.on_off_ratio = self.Id_forward[-1] / np.maximum(self.Id_forward[0], I_noise / W)  # On-off ratio (no units)
        # Off current is minimum(Id, noise floor=1e-12A)
        self.SS = np.diff(self.Vg_forward) / np.diff(
            np.log10(np.absolute(self.Id_forward.astype(float)))) * 1e3  # Subthreshold slope (mV/dec)
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
        gmmax_index = np.where(self.gm==np.amax(self.gm)) #Index of maximum gm
        
        # Calculating Linear Interpolation Vt ##
        gmmax_index = gmmax_index[0]        
        Vg_intercept = self.Vg_forward[gmmax_index] - self.Id_forward[gmmax_index] / np.amax(self.gm)
        self.Vt_lin = Vg_intercept - self.Vds / 2  # Linear extrapolated Vt (V)
        
        # Calculating Max Transconductance Derivative Vt ##
        dgmmax_index = np.where(self.dgm[:-9] == np.amax(self.dgm[:-9]))  # Index of maximum gm
        dgmmax_index, _ = find_peaks(self.dgm, height=0)  # Finds the first local max in dgm
        
        if (dgmmax_index.shape[0] > 0):
            dgmmax_index = dgmmax_index[0]            
            self.Vt_dgm = self.Vg_forward[dgmmax_index]  # Threshold voltage from max derivative of gm (V)
            self.Vov = self.Vg_forward - self.Vt_dgm  # Overdrive voltage wrt dgm Vt (V)
        else:
            self.Vov = self.Vg_forward - self.Vt_lin #Overdrive voltage wrt linear Vt (V)
        # self.Vov = self.Vg_forward - self.Vt_lin #Overdrive voltage wrt linear Vt (V)
        # print(self.Vt_lin)        


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


class IdVd:
    """Class for Id-Vd data"""

    def __init__(self, inp):  # Argument is of type InputData
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

    def plot_IdVd(self, ax, fig, Vgs=None):
        """Plots Id-Vd characteristics in linear scale"""
        flag = 0
        
        if Vgs is None:  # When plot_IdVd is called directly, not through TLM
            flag = 1
            ax.set(xlabel=vd_label, ylabel=id_label)
            ax.annotate(r'$L_{ch}$ = ' + self.L_label, xy=(annotate_x, annotate_y1), xycoords='axes fraction')
            ax.plot(self.Vd, self.Id * 1e6)  # Plotting Id (uA/um) vs Vds (V) for alL Vgs            
        else:            
            i = np.nonzero(np.round(self.Vg, 1) == float(Vgs))                           
            i = i[0].item()  # Finding the index based on Vg
            ax.plot(self.Vd, self.Id[:, i] * 1e6)
            # Plotting Id (uA/um) vs Vds (V) for given Vgs
            
        if flag == 1:            
            ax.locator_params(nbins=6, axis='y')            
            legends = [str(i) + ' V' for i in self.Vg]
            ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.05])
            #Set tick frequency
            set_ticks(ax, 5)
            if isbipolar_idvd == 1: # in case of bipolar Id-Vd, indicate origin
                ax.axhline(0, color='grey', ls='dotted', lw=0.6)
                ax.axvline(0, color='grey', ls='dotted', lw=0.6)                            
        
class TLM:
    """Class for TLM data"""

    def __init__(self, inp_data, inp_count):
        self.data = inp_data  # Input data for all channel lengths
        self.count = inp_count  # Number of channel lengths        
        self.idvg = np.empty(shape=(self.count,num_vd_idvg), dtype=object)
        # self.idvg = np.empty((self.count,num_vg_idvd))  # Id-Vg object
        # self.idvg = [None] * self.count  # Id-Vg object
        self.idvd = np.empty(shape=(self.count,), dtype=object)  # Id-Vd object        
        self.L = np.empty(shape=(self.count,), dtype=float)  # Channel lengths in the TLME
        self.L_label = np.empty(shape=(self.count,), dtype=object)  # Label strings to use for plotting
        
        for i in np.arange(self.count):
            for j in np.arange(num_vd_idvg):                
                self.idvg[i,j] = IdVg(self.data[i], j)  # Id-Vg object
                self.idvg[i,j].idvg_calc()  # Perform calculations on Id-Vg
                # self.idvg[i] = IdVg(self.data[i])  # Id-Vg object
                # self.idvg[i].idvg_calc()  # Perform calculations 
            self.idvd[i] = IdVd(self.data[i])  # Id-Vd object
            self.L[i] = self.idvg[i,0].L
            # self.L[i] = self.idvg[i].L
            
            if self.L[i] < 1:
                self.L_label[i] = str(int(self.L[i]*1000))+r' nm'
            else:
                self.L_label[i] = str(self.L[i])+r' $\mu$m'  
                
        num_vg_forward = int(self.idvg[0,0].num_vg_forward)  # Number of Vg points in forward sweep of Id-Vg data
        # num_vg_forward = int(self.idvg[0].num_vg_forward)  # Number of Vg points in forward sweep of Id-Vg data
        num_vg = int(self.idvg[0,0].num_vg)  # Number of Vg points in Id-Vg data        
        # num_vg = int(self.idvg[0].num_vg)  # Number of Vg points in Id-Vg data        
        self.Vds_idvg = [self.idvg[0,0].Vds,  self.idvg[0,-1].Vds] # Vds forId-Vg sweep        
        self.Vov = np.zeros([num_vg_forward, self.count])  # Overdrive voltage (V)
        self.Vg_idvg = self.idvg[0,0].Vg  # Vg used in Id-Vg (V)        
        self.Id_idvg = np.zeros([num_vg, self.count,num_vd_idvg])        
        self.Vg_idvg_forward = self.idvg[0,0].Vg_forward  # Vg used in Id-Vg (V)                
        # Id from Id-Vg collected for all channel lengths(A/um)
        self.Vg_idvd = self.idvd[0].Vg  # Vg used in Id-Vd (V)
        self.Id_adjusted = []  # Id adjusted after Vt alignment (A/um)
        self.Rtot = []  # Total resistance (Ohm.um)
        self.Rc = []  # Contact resistance (Ohm.um)
        self.del_Rc = []  # Error in Contact resistance (Ohm.um)
        self.Rsh = []  # Sheet resistance (Ohm/sq)
        self.del_Rsh = []  # Error in Sheet resistance (Ohm/sq)
        self.mu_TLM = []  # Mobility extracted from TLM (cm2/V/s)
        self.n2D = []  # 2D carrier concentration (cm-2)
        self.TLM_fit = []  # Pandas dataframe to store regression data @ highest Vov
        self.R_squared = []  # R-squared parameter for TLM fit

    def tlm_calc(self, vd_step = 'low'):
        """Calculates Rtot,Rc,Rsh,mu_TLM
        vd_step = 'low' for low_vds
        vd_step = 'hi' for high_vds"""        
        
        if vd_step == 'low':
            j = 0
        elif vd_step == 'hi':
            j = -1
        # Finding Vgs-Vt overlap and re-arrangingck
        
        for i in np.arange(self.count):
            self.Vov[:, i] = self.idvg[i,j].Vov                                  
            self.Id_idvg[:, i,0] = self.idvg[i,0].Id
            self.Id_idvg[:, i,-1] = self.idvg[i,-1].Id            
            
        self.Vov = np.round(self.Vov)
        Vov_hi = np.amin(self.Vov.max(0))  # Higher and lower limits of Vov (or Vgs-Vt)
        Vov_lo = np.amax(self.Vov.min(0))
        index = np.logical_and(self.Vov <= Vov_hi, self.Vov >= Vov_lo)
        
        # Window of Vgs-Vt overlap across multiple Lchannel's
        self.Vov = self.Vov[np.ix_(index[:, 0], [0])].flatten()
        self.Id_adjusted = np.zeros([self.Vov.shape[0], self.count])
        # Adjusting Id to align the Vgs-Vt for various channels                
        for i in np.arange(self.count):
            self.Id_adjusted[:, i] = self.idvg[i,j].Id_forward[index[:, i]]            
        
        # Calculating Rc, del_Rc, Rsh, del_Rsh
        self.Rtot = self.Vds_idvg[j] / self.Id_adjusted  # Total resistance (Ohm.um)
        self.Rc = np.zeros(self.Vov.shape[0])
        self.del_Rc = np.zeros(self.Vov.shape[0])
        self.Rsh = np.zeros(self.Vov.shape[0])
        self.del_Rsh = np.zeros(self.Vov.shape[0])
        
        # TLM fitting
        if self.count > 1:
            for i in np.arange(self.Vov.shape[0]):
                X = sm.add_constant(self.L.conj().transpose())  # Adding constant term to the model
                model = sm.OLS(self.Rtot[i, :].conj().transpose(), X)
                results = model.fit()
                self.Rc[i] = results.params[0] / 2  # Contact resistance (Ohm.um)
                self.del_Rc[i] = results.bse[0] / 2  # Error in Contact resistance (Ohm.um)
                self.Rsh[i] = results.params[1]  # Sheet resistance (Ohm/sq)
                self.del_Rsh[i] = results.bse[1]  # Error in Sheet resistance (Ohm/sq)
            self.R_squared = results.rsquared
            d = {rtot_label: self.Rtot[i, :].conj().transpose() / 1e3,
                 lch_label: self.L.conj().transpose()}  # Creating dictionary of Rtot vs L
            self.TLM_fit = pd.DataFrame(data=d)  # Converting dictionary to Pandas dataframe        
            self.n2D = Cox * 1e-9 * (self.Vov + np.spacing(1)) / q  # 2D carrier concentration (cm-2)        
            # np.spacing(1) is used to avoid divide-by-zero error in next step
            self.mu_TLM = 1 / ((self.Rsh * self.n2D) * q + np.spacing(1))  # Mobility extracted from TLM (cm^2/V/s)
            self.del_mu_TLM = self.del_Rsh / (self.Rsh + np.spacing(1)) * self.mu_TLM  # Error in mobility extracted from TLM (cm^2/V/s)

    def plot_TLM_fit(self, ax, fig):
        """Plots TLM fit - regression plot"""            
        
        plt.gca().set_prop_cycle(None)
        ax.set_xlim(left=0, right=round(self.L[-1]*1.1/0.2)*0.2)  # Sets lower x-limit to 0
        sns.regplot(x=lch_label, y=rtot_label, data=self.TLM_fit,\
                        ci=95, truncate=False, ax=ax)         
        # truncate = False plots the regression line to fill x-limits
        start, _ = ax.get_ylim()
        ax.set_ylim(bottom=np.minimum(0, start))  # Sets lower y-limit to 0
        Rc_anno = np.round(self.Rc[-1] / 1e3, 1)
        delRc_anno = np.round(self.del_Rc[-1] / 1e3, 1)
        ax.annotate(r'$R_{C}$ = ' + str(Rc_anno) + r' $\pm$ ' + str(delRc_anno) + r' k$\Omega$.$\mu$m',
                    xy=(annotate_x, annotate_y1), xycoords='axes fraction')
        ax.annotate(r'at $n_{2D}$ = ' + str(np.round(self.n2D[-1] / 1e12, 1)) + r' $\times$ $10^{12} cm^{-2}$',
                    xy=(annotate_x, annotate_y22), xycoords='axes fraction', fontsize=8.5)        
        ax.annotate(r'R-squared = ' + str(np.round(self.R_squared, 3)),
                    xy=(.6, .1), xycoords='axes fraction', fontsize=8.5)
        # set_ticks(ax,5)
        set_ticks(ax,5,'x',0.2)

    def plot_mu_TLM(self, ax, fig, flag=None):
        """Plots mu_TLM vs Vov or n2D
            1) flag == None - plot vs Vov
            2) flag == 'n2D'- plot vs n2D"""        
    
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
        ax.annotate(r'at $n_{2D}$ = ' + str(np.round(self.n2D[-1] / 1e12, 1))\
                    + r' $\times$ $10^{12} cm^{-2}$',
                    xy=(annotate_x, annotate_y22), xycoords='axes fraction', fontsize=8.5)
        set_ticks(ax,5)
        
    def plot_Rc(self, ax, fig, flag=None):
        """Plots Rc vs Vov or n2D
            1) flag == None - plot vs Vov
            2) flag == 'n2D'- plot vs n2D"""    
    
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
            ax.set_ylim(bottom=(np.minimum(np.amin(self.Rc[self.n2D >= n2D_LL]\
                                 - self.del_Rc[self.n2D >= n2D_LL]) / 1e3, 0)))            
            # Setting lower x-limit to max(n2D_LL,n2D[0])              
            ax.set_xlim(left=np.maximum(n2D_LL,self.n2D[0])/1e12)
        set_ticks(ax,5)

    def plot_Rsh(self, ax, fig, flag=None):
        """Plots Rsh vs Vov or n2D
            1) flag == None - plot vs Vov
            2) flag == 'n2D'- plot vs n2D"""
    
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
            ax.annotate(r'$V_{DS}$ = ' + str(self.Vds_idvg[j]) + ' V', xy=(.65, .65), xycoords='axes fraction')
            ax.plot(self.Vg_idvg, self.Id_idvg[:,:,j] * 1e6)  # Plotting Id (uA/um) vs Vgs (V)
        elif vd_step == 'hi':            
            j = -1
            ax.annotate(r'$V_{DS}$ = ' + str(self.Vds_idvg[j]) + ' V', xy=(.65, .65), xycoords='axes fraction')
            ax.plot(self.Vg_idvg, self.Id_idvg[:,:,j] * 1e6)  # Plotting Id (uA/um) vs Vgs (V)
        legends = [i for i in self.L_label]
        ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])

    def plot_IdVg_multVds(self, ax, fig, lg='lin',channel=0):
        """
        channel: 0-5 corresponds to CH1 - CH6
        """
        
        self.idvg[channel,0].plot_IdVg(ax, fig, lg)  # Plotting Id-Vg for all channel lengths in log scale
        self.idvg[channel,-1].plot_IdVg(ax, fig, lg)  # Plotting Id-Vg for all channel lengths in log scale
        legends = [r'$V_{DS}$ = ' + str(self.Vds_idvg[0]) + ' V', r'$V_{DS}$ = ' + str(self.Vds_idvg[-1]) + ' V']
        ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])

    def plot_IdVg_Vov(self, ax, fig, lg='lin'):
        """Plots Id-Vg for all channel lengths in log scale
            1)lg == None for linear scale 
            2)lg == 'log' for log scale"""
    
        ax.annotate(r'$V_{DS}$ = ' + str(self.Vds_idvg) + ' V', xy=(.15, annotate_y1), xycoords='axes fraction')
        ax.set(xlabel=vov_label, ylabel=id_label)
        
        if lg == 'log':
            ax.set_yscale('log')
        ax.plot(self.Vov, self.Id_adjusted * 1e6)  # Plotting Id (uA/um) vs Vgs (V)
        legends = [i for i in self.L_label]
        ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])

    def plot_IdVd(self, ax, fig, Vgs):
        """Plots Id-Vd of for all channel lengths at a given Vgs""" 
        
        for i in np.arange(self.count):
            self.idvd[i].plot_IdVd(ax, fig, Vgs)        
        legends = [i for i in self.L_label]
        idmax_anno = int(np.round(np.amax(np.abs(self.idvd[0].Id*1e6))))
        ax.annotate(r'$V_{GS}$ = ' + str(Vgs) + ' V', xy=(0.03, annotate_y1), xycoords='axes fraction')        
        ax.annotate(r'$I_{D,max}$ = ' + str(idmax_anno)\
         + r' $\mu A / \mu m$', xy=(0.03, annotate_y22), xycoords='axes fraction')
        ax.legend(legends, loc='lower right',bbox_to_anchor=[0.95,0.0])
        ax.set(xlabel=vd_label, ylabel=id_label)
        ax.locator_params(nbins=6, axis='y')        
        #Set tick frequency
        set_ticks(ax,5)
            
        if isbipolar_idvd == 1:
            ax.axhline(0, color='grey', ls='dotted', lw=0.6)
            ax.axvline(0, color='grey', ls='dotted', lw=0.6)

# Main section of code
os.chdir(target_dir)
tlm_set = [None] * device_index.shape[0]  # Creating None array to hold TLM objects

for i in device_index:
    if pd.isnull(file_prefix[i]):  # Empty file prefix is read as NaN from excel
        file_prefix[i] = ""
    if isfolder_organized[i]:  # If organized into folders, cd into folder
        os.chdir(device[i])
    count = 0
    col_index = device_params[np.ix_([i], [5, 6, 7, 8, 9, 10])]  # Column indices for various values
    col_index = col_index[0]  # Converting 2-D array to 1-D
    # 1,2,3 - Ids, Vgs, Vds for Id-Vg
    # 4,5,6 - Ids, Vgs, Vds for Id-Vd
    num_vg_idvd = device_params[i, 11]  # Number of Vg points in Id-Vd sweep
    isbipolar_idvd = device_params[i, 12]  # Whether Id-Vd sweep is bipolar
    num_var_idvd = device_params[i, 13]  # Number of variables in a single Id-Vd measurement
    W = device_params[i, 14]  # Channel width (um)
    tox = device_params[i, 15]  # Thickness of the gate oxide (nm)
    input_data = [None] * channel_length.shape[0]  # Creating dummy array for input
    
    for j in channel_index:        
        filename_g = file_prefix[i] + channel_length[j] + "-G.xls"  # Filename of id-vg file
        filename_d = file_prefix[i] + channel_length[j] + "-D.xls"  # Filename of id-vd file
        # Reading Id-Vg files
        my_file_g = Path(filename_g)  # Creating file object        
        my_file_d = Path(filename_d)  # Creating file object
        
        if my_file_g.exists() and my_file_d.exists():  # Checking if IdVg and IdVd files exist            
            wb_idvg = xlrd.open_workbook(filename_g, logfile=open(os.devnull, 'w'))
            id_vg = np.asarray(pd.read_excel(wb_idvg, engine='xlrd'))  # Reading Id-Vg data            
            wb_idvd = xlrd.open_workbook(filename_d, logfile=open(os.devnull, 'w'))
            id_vd = np.asarray(pd.read_excel(wb_idvd, engine='xlrd'))  # Reading Id-Vd data            
            input_data[count] = InputData(id_vg, id_vd, j)  # Assigning the read data to InputData object
            count += 1
    tlm_set[i] = TLM(input_data[:count], count)  # Creating TLM object
    # tlm_set[i].tlm_calc(vd_step = 'low')  # Performing TLM calculations    
    tlm_set[i].tlm_calc(vd_step = 'hi')  # Performing TLM calculations    
    
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
    tlm_set[i].plot_IdVg_multVds(ax3, fig, 'log', channel = 0)  # Plotting Id-Vg for 100nm in log scale for smallest and largest vds
    # tlm_set[i].plot_IdVd(ax2, fig, Vg_idvd_plt)  # Plotting Id-Vd for all channel lengths @ Vgs=50V
    # tlm_set[i].idvd[0].plot_IdVd(ax3, fig)  # Plotting Id-Vd for 1st channel length
    
    if tlm_set[i].count > 1:
        tlm_set[i].plot_TLM_fit(ax4, fig)  # Plotting TLM fit for the last Vov
        tlm_set[i].plot_Rc(ax5, fig, flag = 'n2D')  # Plotting Rc vs n2D
        tlm_set[i].plot_mu_TLM(ax6, fig, flag = 'n2D')  # Plotting mu_eff vs n2D
        # tlm_set[i].plot_Rsh(ax7, fig, 'n2D')  # Plotting Rsh vs n2D
        rsh = tlm_set[i].Rsh[-1]
        rc = tlm_set[i].Rc[-1]
        
    plt.show()
    fig.savefig("summary-" + str(i) + ".svg", transparent=True,\
                                        bbox_inches='tight', pad_inches=0.1)    
    if isfolder_organized[i]:
        os.chdir(target_dir)