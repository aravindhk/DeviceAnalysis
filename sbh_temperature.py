# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 12:23:25 2020

@author: Aravindh Kumar
"""
"""
Code to analyse temperature dependent data from Janis
"""

#import sys
import os
#Including 3rd-party libraries
#import matplotlib
#matplotlib.rcParams['text.usetex'] = True #For rendering latex in plot elements
import matplotlib.pyplot as plt
import matplotlib as mpl
#from scipy.signal import find_peaks  #For finding local maxima in sbh
#Changing plot properties
plt.style.use(['default', 'seaborn-bright'])
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['xtick.labelsize'] = 9
mpl.rcParams['ytick.labelsize'] = 9
mpl.rcParams['font.size'] = 10
mpl.rcParams['lines.markersize'] = 3
mpl.rcParams['figure.subplot.wspace'] = 0.3 #Horizontal spacing in subplots
mpl.rcParams['figure.subplot.hspace'] = 0.3 #Vertical spacing in subplots
mpl.rcParams['figure.figsize'] = 12,3 #Figure size
mpl.rcParams['figure.constrained_layout.h_pad'] =  0.04167  #Padding around axes objects
mpl.rcParams['figure.constrained_layout.w_pad'] =  0.04167  
mpl.rcParams['legend.fontsize'] = 'small'
#mpl.rcParams['figure.autolayout'] = True

#plt.style.use('default')
import pandas as pd
import xlrd
#import csv
import numpy as np
import seaborn as sns
#import scipy.linalg
import statsmodels.api as sm
#from statsmodels.api import graphics
from pathlib import Path

#Setting current directory and target directory
current_dir = os.getcwd()
user_folder = os.getenv('USERPROFILE');
##### Selecting the folder ####################################################
target_dir = user_folder + "\\Google Drive\\Research\\Projects\\PSG Doping\\Janis\\2020-03-10-P12-02-Temperature"
#target_dir = user_folder + "\\Google Drive\\Research\\Projects\\PSG Doping\\Janis\\2020-03-10-P13-01-Temperature"
###############################################################################

os.chdir(target_dir) #Move to target directory

#Reading file parameters from input files
#Reading channel parameters - how many channels, what channel names, etc.
temp_params = Path("temperature_parameters.xlsx") #Channel parameters file
temp_params = xlrd.open_workbook(temp_params,logfile=open(os.devnull,'w'))
temp_params = np.asarray(pd.read_excel(temp_params,engine='xlrd')) #Reading Id-Vg data
temperature = temp_params[:,0] #Channel lengths (um)
temperature_index = np.arange(temperature.shape[0]) #Channel lengths index

#Reading device parameters
device_params = Path("device_parameters.xlsx") #Device parameters file - what sweeps, column indices, etc.
device_params = xlrd.open_workbook(device_params,logfile=open(os.devnull,'w'))
device_params = np.asarray(pd.read_excel(device_params,engine='xlrd')) #Reading Id-Vg data
chip = device_params[:,2] #List of chips (needed?)
device = device_params[:,2] #List of devices
file_prefix = device_params[:,3] #List of file prefixes
isfolder_organized = device_params[:,4] #Whether the device data is organized into folders
col_index = device_params[np.ix_([0],[5,6,7,8,9,10])] #Column indices for various values
col_index = col_index[0] #Converting 2-D array to 1-D
##1,2,3 - Ids, Vgs, Vds for Id-Vg
##4,5,6 - Ids, Vgs, Vds for Id-Vd
#num_vg_idvd = device_params[0,11] #Number of Vg points in Id-Vd sweep
#isbipolar_idvd = device_params[0,12] #Whether Id-Vd sweep is bipolar
#num_var_idvd = device_params[0,13] #Number of variables in a single Id-Vd measurement
W = device_params[0,14] #Channel width (um)
tox = device_params[0,15] #Thickness of the gate oxide (nm)
q = 1.6e-19 #Charge of an electron (C)
Cox_30nm = 116 #30nm SiO2 capacitance (Kirby ACS Nano) (nF/cm2)
Cox = Cox_30nm*30/tox #100nm SiO2 capacitance (nF/cm2)
Id_thresh = 2e-8 #For constant Vt calculation (A)
I_noise = 1e-12 #Noise floor for Id (A)
n2D_LL = 2e12 #n2D lower limit for plotting (cm-2)
kB1000 = 86.7 #kB*T at T = 1000 K (meV)

chip_index = np.arange(chip.shape[0]) #Chip index
device_index = np.arange(device.shape[0]) #Device index


def smooth(a,WSZ=7):
#     a: NumPy 1-D array containing the data to be smoothed
#     WSZ: smoothing window size needs, which must be odd number,
#     as in the original MATLAB implementation
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ    
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , out0, stop  ))
###############################################################################

class InputData:
    """Class for input data"""
    def __init__(self,idvg,index):
        self.id_vg= idvg #Assigning read Id-Vg data        
        self.temperature = temperature[index] #Channel length (value)
#        self.id_vd = idvd #Assigning read Id-Vd data
###############################################################################

class IdVg:
    """Class for Id-Vg data"""
    def __init__(self,inp): #Argument is of type InputData        
        self.T = inp.temperature #Temperatures (K)        
        self.inv_T = 1000/inp.temperature #Inverse of temperatures (K-1)
        self.num_vg = int(inp.id_vg.shape[0]) #Number of Vg points in Id-Vg
        self.num_vg_forward = int(inp.id_vg.shape[0]/2) #Number of Vg points in 
#                                                forward sweep for calculations       
#        self.Id = inp.id_vg[0:self.num_vg,8]/W #Parsing only forward sweep drain current (A/um)
#        self.Vg = inp.id_vg[0:self.num_vg,4] #Parsing only forward sweep gate voltage (V)        
#        self.Vds = inp.id_vg[0,9] #Reading Vds from Id-Vg data (V)
        self.Id = inp.id_vg[0:self.num_vg,col_index[0]]/W #Reading drain current (A/um)
        self.Vg = inp.id_vg[0:self.num_vg,col_index[1]] #Reading gate voltage (V)      
        self.Id_forward  = inp.id_vg[0:self.num_vg_forward,col_index[0]]/W 
        #Parsing only forward sweep drain current (A/um)
        self.Vg_forward = inp.id_vg[0:self.num_vg_forward,col_index[1]] 
        #Parsing only forward sweep gate voltage (V)      
        self.Vds = np.round(inp.id_vg[0,col_index[2]],1) #Reading Vds from Id-Vg data (V)
#        self.Vds = 1 #Hard-coding Vds for Id-Vg data (V)
        self.I_off = np.maximum(self.Id[0],I_noise) #Finding off current 
                                                #(comparing with noise floor)        
        self.on_off_ratio = [] #On-off ratio (no units)
        self.SS = [] #Subthreshold swing (mV/dec)
#        self.mu_FE = [] #Field-effect mobility (cm2/V/s)
        self.gm = [] #Transconductance (S/um)
        self.dgm = [] #Transconductance derivative (S/V/um)
        self.Vt_const = [] #Constant current threshold voltage (V)
        self.Vt_lin = [] #Linear extrapolation threshold voltage (V)
        self.Vt_dgm = [] #Threshold voltage from max derivative of gm (V)
        self.Vov = [] #Overdrive voltage Vgs-Vt (V)
        
    def idvg_calc(self):
        """Calculation of parameters from Id-Vg"""
        self.on_off_ratio = self.Id_forward[-1]/np.maximum(self.Id_forward[0],I_noise/W) #On-off ratio (no units)
#        (Off current is determined by comparing minimum Id with noise floor which is 1e-12 A)                
        self.SS = np.diff(self.Vg_forward)/(np.diff(np.log10(np.absolute(self.Id_forward.astype(float))))*1e3+np.spacing(1)) #Subthreshold slope (mV/dec)
        #Calculating field-effect mobility                
        self.mu_FE = np.zeros((self.Vg_forward.shape[0]-1))      
        self.gm = np.diff(smooth(self.Id_forward,5))/np.diff(self.Vg_forward) #Transconductance (S/um)
        self.dgm = np.diff(smooth(self.gm,5))/np.diff(self.Vg_forward[:-1]) #Transconductance derivative (S/V/um)
        self.dgm = smooth(self.dgm,5) #Smoothed Transconductance derivative (S/V/um)
        self.gm = smooth(self.gm) #Smoothed Transconductance (S/um)
#        self.mu_FE = self.L*self.gm/(Cox*self.Vds*1e-9) #Possible under-estimation? (cm2/V/s)
#         Calculating Constant Current Vt ##
#        if np.amin(np.absolute(self.Id_forward)) < Id_thresh*W/self.L:
#            self.Vt_const = np.amax(self.Vg_forward[self.Id_forward<Id_thresh*W/self.L]) #Constant current Vt (V)
        gmmax_index = np.where(self.gm==np.amax(self.gm)) #Index of maximum gm
#         Calculating Linear Interpolation Vt ##
        Vg_intercept = self.Vg_forward[gmmax_index] - self.Id_forward[gmmax_index]/np.amax(self.gm)
        self.Vt_lin = Vg_intercept - self.Vds/2 #Linear extrapolated Vt (V)
#         Calculating Max Transconductance Derivative Vt ##
        dgmmax_index = np.where(self.dgm[:-9]==np.amax(self.dgm[:-9])) #Index of maximum gm
        self.Vt_dgm = self.Vg_forward[dgmmax_index] #Threshold voltage from max derivative of gm (V)    
#        self.Vov = self.Vg - self.Vt_lin #Overdrive voltage wrt linear Vt (V)        
        self.Vov = self.Vg_forward - self.Vt_dgm #Overdrive voltage wrt dgm Vt (V) 
#        print(self.Vt_dgm)
        
    def plot_IdVg(self,ax,fig=None,lg='lin'):
        """Plots Id-Vg characteristics. 
        1) lg == 'lin' for linear scale 
        2) lg == 'log' for log scale"""                    
#        fig, ax = plt.subplots()
        ax.plot(self.Vg, self.Id*1e6) #Plotting Id (uA/um) vs Vgs (V)                
        if fig != None: #This happens when plot_IdVg is called directly, not through TLM
            ax.annotate(r'$V_{DS}$ = '+str(self.Vds)+' V',xy=(.45, .85), xycoords='figure fraction')
            ax.set(xlabel='Gate Voltage (V)', ylabel=r'Drain Current ($\mu A$/$\mu m$)')
#            ax.grid()                
            if lg == 'lin':
                ax.locator_params(nbins=6, axis='y')
#            plt.ylim(-5, 20*np.ceil(1e6*np.amax(self.Id)/20))
#            Ensures major tick is at plot edge
            elif lg == 'log':
                ax.set_yscale('log') 
#            plt.ylim(10**np.floor(np.log10(1e6*np.amin(np.absolute(self.Id)))),
#                     10**np.ceil(np.log10(1e6*np.amax(np.absolute(self.Id))))) 
#            Ensures major tick is at plot edge
#        ax.plot(self.Vg, self.Id*1e6) #Plotting Id (uA/um) vs Vgs (V)        
###############################################################################
#
#class IdVd:
#    """Class for Id-Vd data"""
#    def __init__(self,inp): #Argument is of type InputData
#        if isbipolar_idvd: #If bipolar sweep, just take forward sweep, else dual sweep
#            num_vd = int(inp.id_vd.shape[0])
#        else:
#            num_vd = int(inp.id_vd.shape[0])
#        self.Vd = np.round(inp.id_vd[0:num_vd,col_index[5]],2) #Drain voltage (V)
#        self.Id = inp.id_vd[np.arange(num_vd),col_index[3]:num_vg_idvd*num_var_idvd:num_var_idvd]/W #Drain current (A/um)
#        self.Vg = inp.id_vd[0,col_index[4]:num_vg_idvd*num_var_idvd:num_var_idvd]
#        self.temperatureength = inp.temperatureength #Channel length string
#        self.L = inp.temperature #Channel length (um)
#        self.Id_max = np.amax(np.abs(inp.id_vd[:,col_index[3]+num_var_idvd*(num_vg_idvd-1)]/W*1e6)) #Maximum drain current (uA/um)
#    
#    def plot_IdVd(self,ax,fig,Vgs=None):
#        """Plots Id-Vd characteristics in linear scale"""
#        flag = 0
#        if Vgs is None: #This happens when plot_IdVd is called directly, not through TLM
#            flag = 1
#            fig, ax = plt.subplots()                                 
#            ax.set(xlabel='Drain Voltage (V)', ylabel=r'Drain Current ($\mu A$/$\mu m$)') 
#            ax.annotate(r'$L_{CH}$ = '+str(self.L)+r' $\mu m$',xy=(.45, .85), xycoords='axes fraction')
#            ax.plot(self.Vd, self.Id*1e6) #Plotting Id (uA/um) vs Vds (V) for alL Vgs                  
#            if isbipolar_idvd == 1:
#                ax.axhline(0, color='grey',ls='dotted',lw=0.6)
#                ax.axvline(0, color='grey',ls='dotted',lw=0.6)
#        else:            
#            i = np.nonzero(self.Vg==float(Vgs))
#            i = i[0].item() #Finding the index based on Vg
#            ax.plot(self.Vd, self.Id[:,i]*1e6)             
#            #Plotting Id (uA/um) vs Vds (V) for given Vgs                                
#        if flag == 1:                    
#            for k in np.arange(self.Vg.shape[0]):                
#                text_offset = (-50,3)
#                if k == self.Vg.shape[0]-1:
#                    text_offset = (-50,-35)
#                plt.annotate(str(self.Vg[k])+' V', 
#                             xy=(np.amax(self.Vd), np.amax(self.Id[:,k]*1e6)), xycoords='data',
#                             xytext=text_offset, textcoords='offset points', 
#                             fontsize=16)
#            ax.locator_params(nbins=6, axis='y')
#            plt.ylim(-5, 50*np.ceil(1e6*np.amax(self.Id)/50)) #Ensures major tick is at plot edge
#            legends = [str(i)+' V' for i in self.Vg]                        
#            ax.legend(legends,loc='best') 
###############################################################################
            
class TLM:
    """Class for TLM data"""    
    def __init__(self,inp_data,inp_count):        
        self.data = inp_data #Input data for all channel lengths              
        self.count= inp_count #Number of channel lengths
        self.idvg = [None]*self.count #Id-Vg object
#        self.idvd = [None]*self.count #Id-Vd object
        self.T = np.zeros(self.count) #Temperatures in the measurement (K)
        self.inv_T = np.zeros(self.count) #Inverse T in the measurement (K-1)
        for i in np.arange(self.count):            
            self.idvg[i] = IdVg(self.data[i]) #Id-Vg object            
            self.idvg[i].idvg_calc() #Perform calculations on Id-Vg
#            self.idvd[i] = IdVd(self.data[i]) #Id-Vd object                
            self.T[i] = self.idvg[i].T
            self.inv_T[i] = self.idvg[i].inv_T
        num_vg_forward = int(self.idvg[0].num_vg_forward) #Number of Vg points in forward sweep of Id-Vg data
#        num_vg= int(self.idvg[0].num_vg) #Number of Vg points in Id-Vg data
        self.Vds_idvg = self.idvg[0].Vds #Vds forId-Vg sweep
        self.Vov = np.zeros([num_vg_forward,self.count]) #Overdrive voltage (V)
#        self.Vg_idvg = self.idvg[0].Vg #Vg used in Id-Vg (V)
#        self.Id_idvg = np.zeros([num_vg,self.count]) 
#        self.Vg_idvg_forward = self.idvg[0].Vg_forward #Vg used in Id-Vg (V)
#        self.Id_idvg_forward = np.zeros([num_vg_forward,self.count]) 
#        Id from Id-Vg collected for all channel lengths(A/um)        
#        self.Vg_idvd = self.idvd[0].Vg #Vg used in Id-Vd (V)
        self.Id_adjusted = [] #Id adjusted after Vt alignment (A/um)
        self.Rtot = [] #Total resistance (Ohm.um)
#        self.Rc = [] #Contact resistance (Ohm.um)
#        self.del_Rc = [] #Error in Contact resistance (Ohm.um)
        self.sbh = [] #Sheet resistance (Ohm/sq)
        self.del_sbh = [] #Error in Sheet resistance (Ohm/sq)
#        self.mu_TLM = [] #Mobility extracted from TLM (cm2/V/s)
        self.n2D = [] #2D carrier concentration (cm-2)
        self.TLM_fit = [] #Pandas dataframe to store regression data @ highest Vov
        self.R_squared = [] #R-squared parameter for TLM fit

    def tlm_calc(self):     
        """Calculates Rtot,Rc,sbh,mu_TLM"""
#        Finding Vgs-Vt overlap and re-arranging
        Vov_max = np.zeros(self.count)
        Vov_min = np.zeros(self.count)
        index = [None]*self.count
        for i in np.arange(self.count):
            Vov_max[i] = self.idvg[i].Vov.max(0)
            Vov_min[i] = self.idvg[i].Vov.min(0)
#            self.Vov[:,i] = self.idvg[i].Vov 
#             self.Id_idvg_forward[:,i] = self.idvg[i].Id_forward
#            self.Id_idvg[:,i] = self.idvg[i].Id
#        self.Vov = np.round(self.Vov)
        Vov_hi = np.amin(Vov_max) #Higher and lower limits of Vov (or Vgs-Vt)
        Vov_lo = np.amax(Vov_min)
        for i in np.arange(self.count):
            index[i] = np.logical_and(np.round(self.idvg[i].Vov)<Vov_hi,np.round(self.idvg[i].Vov)>Vov_lo)
#                            Window of Vgs-Vt overlap across multiple Lchannel's
        self.Vov = np.round(self.idvg[0].Vov[index[0]].flatten())
        self.Id_adjusted = np.zeros([self.Vov.shape[0],self.count])
#                            Adjusting Id to align the Vgs-Vt for various channels
        for i in np.arange(self.count):
             self.Id_adjusted[:,i] = self.idvg[i].Id_forward[index[i]]             
        self.Rtot = self.Vds_idvg/self.Id_adjusted #Total resistance (Ohm.um)
#        Calculating Rc, del_Rc, sbh, del_sbh
#        self.Rc = np.zeros(self.Vov.shape[0])
#        self.del_Rc = np.zeros(self.Vov.shape[0])
        self.sbh = np.zeros(self.Vov.shape[0])
        self.del_sbh = np.zeros(self.Vov.shape[0])
        self.R_squared = np.zeros(self.Vov.shape[0])
        self.TLM_fit = [None]*self.Vov.shape[0]
##### TLM fitting #############################################################
        for i in np.arange(self.Vov.shape[0]):
            X = sm.add_constant(self.inv_T.conj().transpose()) #Adding constant term to the model
            log_Ids = np.log(np.abs(self.Id_adjusted[i,:].conj().transpose()))
            model = sm.OLS(log_Ids,X)         
            results = model.fit()            
#            self.Rc[i] = results.params[0]/2 #Contact resistance (Ohm.um)
#            self.del_Rc[i] = results.bse[0]/2 #Error in Contact resistance (Ohm.um)
            self.sbh[i] = -kB1000*results.params[1] #SBH (meV)            
            self.del_sbh[i] = kB1000*results.bse[1] #Error in SBH (meV)            
#        print(results.summary())
#        print(results.rsquared)
            self.R_squared[i] = results.rsquared
            d = {r'log($I_{DS}$)': log_Ids, 
             r'1000/T ($K^{-1}$)': self.inv_T.conj().transpose()} #Creating dictionary of Rtot vs L
            self.TLM_fit[i] = pd.DataFrame(data=d) #Converting dictionary to Pandas dataframe        
        self.sbh = smooth(self.sbh,5) #Smoothed SBH (meV)
        self.del_sbh = smooth(self.del_sbh,5) #Smoothed Error SBH (meV)
        self.d_sbh = np.diff(self.sbh)/np.diff(self.Vov) #Derivative of SBH
        self.d_sbh = smooth(self.d_sbh,5) #Smoothed Derivative of SBH
        self.d_d_sbh = np.diff(self.d_sbh)/np.diff(self.Vov[:-1]) #2nd Derivative of SBH
        self.d_d_sbh = smooth(self.d_d_sbh,5) #Smoothed 2nd Derivative of SBH            
###############################################################################
        self.n2D = Cox*1e-9*(self.Vov+np.spacing(1))/q #2D carrier concentration (cm-2)
        #np.spacing(1) is used to avoid divide-by-zero error in next step
#        self.mu_TLM = 1/((self.sbh*self.n2D)*q+np.spacing(1)) #Mobility extracted from TLM (cm^2/V/s)
#        self.del_mu_TLM = self.del_sbh/self.sbh*self.mu_TLM #Error in mobility extracted from TLM (cm^2/V/s)
                
    def plot_TLM_fit(self,ax,fig,i):#Vov):
        """Plots TLM fit @ Vov - regression plot"""
#        i = np.nonzero(self.Vov==float(Vov))
#        i = i[0].item() #Finding the index based on Vg
        sns.regplot(x=r'1000/T ($K^{-1}$)', y=r'log($I_{DS}$)', data=self.TLM_fit[i], ci=95, ax=ax)
#        ax.set_ylim(bottom=0) #Sets lower y-limit to 0            
        
#    def plot_mu_TLM(self,ax,fig,flag=None):
#        """Plots mu_TLM vs Vov or n2D
#            1) flag == None - plot vs Vov
#            2) flag == 'n2D'- plot vs n2D"""
##        fig, ax = plt.subplots()        
##        ax.grid()                
#        if flag is None:            
#            ax.set(xlabel='Overdrive Voltage (V)', 
#                   ylabel=r'Contact Resistance (k$\Omega.\mu m$)')#, title=r'$R_C$ vs $V_{GS}-V_{T}$')
#            ax.errorbar(self.Vov[self.Vov>=0], self.mu_TLM[self.Vov>=0],
#                        self.del_mu_TLM[self.Vov>=0]) 
##            Plotting mu_TLM (cm2/V/s) vs Vov(V)
#        elif flag == 'n2D':
#            ax.set(xlabel=r'Carrier concentration ($10^{12}$ $cm^{-2}$)', 
#                   ylabel=r'Effective Mobility ($cm^2$/V/s)')#, title=r'$R_C$ vs $n_{2D}$')
#            ax.errorbar(self.n2D[self.n2D>=n2D_LL]/1e12, self.mu_TLM[self.n2D>=n2D_LL],
#                        self.del_mu_TLM[self.n2D>=n2D_LL]) 
#        ax.set_ylim(bottom=0) #Sets lower y-limit to 0            
##        plt.ylim(5*np.floor(np.amin(self.Rc[self.Vov>=0]/1e3)/5), 
##                 10*np.ceil(np.amax(self.Rc[self.Vov>=0]/1e3)/10))
##            Plotting mu_TLM (cm2/V/s) vs n2D(cm-2)
        
#    def plot_Rc(self,ax,fig,flag=None):
#        """Plots Rc vs Vov or n2D
#            1) flag == None - plot vs Vov
#            2) flag == 'n2D'- plot vs n2D"""
##        fig, ax = plt.subplots()        
##        ax.grid()                
#        if flag is None:            
#            ax.set(xlabel='Overdrive Voltage (V)', 
#                   ylabel=r'Contact Resistance (k$\Omega.\mu m$)')#, title=r'$R_C$ vs $V_{GS}-V_{T}$')
##            ax.errorbar(self.Vov[self.Vov>=0], self.Rc[self.Vov>=0]/1e3,
##                        self.del_Rc[self.Vov>=0]/1e3) 
#            #Plotting Rc (kOhm.um) vs Vov(V)
#            ax.plot(self.n2D[self.n2D>=n2D_LL]/1e12, self.Rc[self.n2D>=n2D_LL]/1e3)             
#        elif flag == 'n2D':
#            ax.set(xlabel=r'Carrier concentration ($10^{12}$ $cm^{-2}$)', 
#                   ylabel=r'Contact Resistance (k$\Omega.\mu m$)')#, title=r'$R_C$ vs $n_{2D}$')
##            ax.errorbar(self.n2D[self.n2D>=n2D_LL]/1e12, self.Rc[self.n2D>=n2D_LL]/1e3,
##                        self.del_Rc[self.n2D>=n2D_LL]/1e3)             
#            ax.plot(self.n2D[self.n2D>=n2D_LL]/1e12, self.Rc[self.n2D>=n2D_LL]/1e3)             
##            ax.set_ylim(bottom=(np.minimum(np.amin(self.Rc-self.del_Rc)/1e3,0))) 
##                                                Sets lower y-limit to 0            
##        plt.ylim(5*np.floor(np.amin(self.Rc[self.Vov>=0]/1e3)/5), 
##                 10*np.ceil(np.amax(self.Rc[self.Vov>=0]/1e3)/10))
##            Plotting Rc (kOhm.um) vs n2D(cm-2)
        
    def plot_sbh(self,ax,fig,order=0,flag=None):        
        """Plots sbh vs Vov or n2D
            1) flag == None - plot vs Vov
            2) flag == 'n2D'- plot vs n2D"""
#        fig, ax = plt.subplots()
#        ax.grid()                
        sbh_max_index = np.where(self.sbh==np.amax(self.sbh)) #Index of maximum gm                        
#        sbh_max_index = 0#sbh_max_index[0][0]
        sbh_max_index = sbh_max_index[0][0]
        if order == 0:
            X_Vov = self.Vov[sbh_max_index:]
            X_n2D = self.n2D[sbh_max_index:]
            Y = self.sbh[sbh_max_index:]
            error = self.del_sbh[sbh_max_index:]
#            error = Y*0
        elif order == 1:
            X_Vov = self.Vov[sbh_max_index:-1]
            X_n2D = self.n2D[sbh_max_index:-1]
            Y = self.d_sbh[sbh_max_index:]
            error = Y*0
        elif order == 2:
            X_Vov = self.Vov[sbh_max_index:-2]
            X_n2D = self.n2D[sbh_max_index:-2]
            Y = self.d_d_sbh[sbh_max_index:]
            error = Y*0
        if flag is None:            
            ax.set(xlabel='Overdrive Voltage (V)', ylabel=
                   r'Schottky barrier height (meV)')#, title=r'$R_{SH}$ vs $V_{GS}-V_{T}$')            
            ax.errorbar(X_Vov, Y, error) 
#            Plotting sbh (kOhm/sq) vs Vov(V)
        elif flag == 'n2D':            
            ax.set(xlabel=r'Carrier concentration ($10^{12}$ $cm^{-2}$)', 
                   ylabel=r'Schottky barrier height (meV)')#, title=r'$R_{SH}$ vs $n_{2D}$')
            ax.errorbar(X_n2D/1e12, Y,error) 
#            Plotting sbh (kOhm/sq) vs n2D(cm-2)
#        ax.locator_params(nbins=6, axis='y')
        ax.set_ylim(bottom=(np.minimum(np.amin(Y-error),0))) 
#                                                Sets lower y-limit to 0          
#        plt.ylim(25*np.floor(np.amin(self.sbh[self.Vov>=0]/1e3)/25), 
#                 25*np.ceil(np.amax(self.sbh[self.Vov>=0]/1e3)/25))
#        Ensures major tick is at plot edge
        
#    def plot_Rtot(self,ax,fig,flag=None):
#        """Plots Rtot vs Vov or n2D
#            1) flag == None - plot vs Vov
#            2) flag == 'n2D'- plot vs n2D"""
##        fig, ax = plt.subplots()
##        ax.grid()                
#        if flag is None:            
#            ax.set(xlabel='Overdrive Voltage (V)', ylabel=
#                   r'Total Resistance (k$\Omega.\mu $m)')#, title=r'$R_{SH}$ vs $V_{GS}-V_{T}$')
#            ax.plot(self.Vov[self.Vov>=0], self.Rtot[self.Vov>=0]/1e3) 
##            Plotting sbh (kOhm/sq) vs Vov(V)
#        elif flag == 'n2D':
#            ax.set(xlabel=r'Carrier concentration ($10^{12}$ $cm^{-2}$)', 
#                   ylabel=r'Total Resistance (k$\Omega.\mu $m)')#, title=r'$R_{SH}$ vs $n_{2D}$')
#            ax.plot(self.n2D[self.n2D>=n2D_LL]/1e12, self.Rtot[self.n2D>=n2D_LL]/1e3) 
##            Plotting sbh (kOhm/sq) vs n2D(cm-2)
#        ax.locator_params(nbins=6, axis='y')
#        ax.set_ylim(bottom=0) #Sets lower y-limit to 0            
##        plt.ylim(25*np.floor(np.amin(self.sbh[self.Vov>=0]/1e3)/25), 
##                 25*np.ceil(np.amax(self.sbh[self.Vov>=0]/1e3)/25))
##        Ensures major tick is at plot edge
        
    def plot_Rc_Rsh(self,ax,fig):
        """Plots Rc and Rsh in the same figure"""        
        
    def plot_IdVg(self,ax,fig,lg='lin'):
        """Plots Id-Vg for all channel lengths in log scale
            1)lg == None for linear scale 
            2)lg == 'log' for log scale"""        
#        fig, ax = plt.subplots()
        for i in np.arange(self.count):
            self.idvg[i].plot_IdVg(ax)
        ax.annotate(r'$V_{DS}$ = '+str(self.Vds_idvg)+' V',xy=(.15, .85), xycoords='axes fraction')
        ax.set(xlabel='Gate Voltage (V)', ylabel=r'Drain Current ($\mu A$/$\mu m$)')
#            ax.grid()                
#        if lg == 'lin':
#            ax.locator_params(nbins=6, axis='y')
#            plt.ylim(-5, 20*np.ceil(1e6*np.amax(self.Id_idvg)/20))
            #Ensures major tick is at plot edge
        if lg == 'log':
            ax.set_yscale('log') 
#            plt.ylim(10**np.floor(np.log10(1e6*np.amin(np.absolute(self.Id_idvg)))),
#                     10**np.ceil(np.log10(1e6*np.amax(np.absolute(self.Id_idvg))))) 
#            Ensures major tick is at plot edge        
#        ax.plot(self.Vg_idvg, self.Id_idvg*1e6) #Plotting Id (uA/um) vs Vgs (V)        
#        legends = [r'$L_{CH}$ = '+str(i)+r' $\mu$m' for i in self.L]                        
        legends = [str(int(i))+r' K' for i in self.T]                                
        ax.legend(legends,loc='best')        

    def plot_IdVg_Vov(self,ax,fig,lg='lin'):
        """Plots Id-Vg for all channel lengths in log scale
            1)lg == None for linear scale 
            2)lg == 'log' for log scale"""        
#        fig, ax = plt.subplots()
        ax.annotate(r'$V_{DS}$ = '+str(self.Vds_idvg)+' V',xy=(.35, .35), xycoords='axes fraction')
        ax.set(xlabel='Overdrive Voltage (V)', ylabel=r'Drain Current ($\mu A$/$\mu m$)')
#            ax.grid()                
#        if lg == 'lin':
#            ax.locator_params(nbins=6, axis='y')
#            plt.ylim(-5, 20*np.ceil(1e6*np.amax(self.Id_idvg)/20))
            #Ensures major tick is at plot edge
        if lg == 'log':
            ax.set_yscale('log') 
#            plt.ylim(10**np.floor(np.log10(1e6*np.amin(np.absolute(self.Id_idvg)))),
#                     10**np.ceil(np.log10(1e6*np.amax(np.absolute(self.Id_idvg))))) 
#            Ensures major tick is at plot edge        
        ax.plot(self.Vov, self.Id_adjusted*1e6) #Plotting Id (uA/um) vs Vgs (V)        
#        legends = [r'$L_{CH}$ = '+str(i)+r' $\mu$m' for i in self.L]                        
        legends = [str(int(i))+r' K' for i in self.T]                                
        ax.legend(legends,loc='best')        
        
#    def plot_IdVd(self,ax,fig,Vgs):
#        """Plots Id-Vd of for all channel lengths at a given Vgs"""
##        fig, ax = plt.subplots()              
#        for i in np.arange(self.count):
#            self.idvd[i].plot_IdVd(ax,fig,Vgs)
#        legends = [r'$L_{CH}$ = '+str(i)+r' $\mu$m' for i in self.L]                
#        legends = [str(i)+r' $\mu$m' for i in self.L]                
#        ax.annotate(r'$V_{GS}$ = '+str(Vgs)+' V',xy=(.45, .85), xycoords='axes fraction')
#        ax.legend(legends,loc='best')        
#        ax.set(xlabel='Drain Voltage (V)', ylabel=r'Drain Current ($\mu A$/$\mu m$)')#, title='Output Characteristics')                
#        ax.grid()
#        plt.ylim(-5, 25*np.ceil(1e6*np.amax(
#                self.idvd[0].Id[:,np.where(self.Vg_idvd==Vgs)[0][0]]/25))) 
#        #Ensures major tick is at plot edge
#        ax.locator_params(nbins=6, axis='y')        
#        if isbipolar_idvd == 1:
#                ax.axhline(0, color='grey',ls='dotted',lw=0.6)
#                ax.axvline(0, color='grey',ls='dotted',lw=0.6)
         
#### Main section of code #####################################################
os.chdir(target_dir)
tlm_set = [None]*device_index.shape[0] #Creating None array to hold TLM objects
for i in device_index:    
    if pd.isnull(file_prefix[i]): #Empty file prefix is read as NaN from excel
        file_prefix[i] = ""
    if isfolder_organized[i]: #If organized into folders, cd into folder
        os.chdir(device[i])    
    count = 0
    col_index = device_params[np.ix_([i],[5,6,7,8,9,10])] #Column indices for various values
    col_index = col_index[0] #Converting 2-D array to 1-D
#    1,2,3 - Ids, Vgs, Vds for Id-Vg
#    4,5,6 - Ids, Vgs, Vds for Id-Vd
#    num_vg_idvd = device_params[i,11] #Number of Vg points in Id-Vd sweep
#    isbipolar_idvd = device_params[i,12] #Whether Id-Vd sweep is bipolar
#    num_var_idvd = device_params[i,13] #Number of variables in a single Id-Vd measurement
    W = device_params[i,14] #Channel width (um)
    tox = device_params[i,15] #Thickness of the gate oxide (nm)
    input_data = [None]*temperature.shape[0] #Creating dummy array for input
    for j in temperature_index:            
        filename_g = file_prefix[i] + str(temperature[j]) + "K-G.xls" #Filename of id-vg file        
#        filename_d = file_prefix[i] + temperature[j] + "-D.xls" #Filename of id-vd file
#        Reading Id-Vg files
        my_file_g = Path(filename_g) #Creating file object
#        my_file_d = Path(filename_d) #Creating file object
        if my_file_g.exists(): #Checking if IdVg and IdVd files exist
            wb_idvg = xlrd.open_workbook(filename_g,logfile=open(os.devnull,'w'))
            id_vg = np.asarray(pd.read_excel(wb_idvg,engine='xlrd')) #Reading Id-
            input_data[count] = InputData(id_vg,j) #Assigning the read data to InputData object
            count += 1        
    tlm_set[i] = TLM(input_data[:count],count) #Creating TLM object
    tlm_set[i].tlm_calc() #Performing TLM calculations
    
#### Plot options #############################################################
#    tlm_set[i].idvg[0].plot_IdVg('log') #Plotting Id-Vg for 1st channel length
#    tlm_set[i].idvd[0].plot_IdVd() #Plotting Id-Vd for 1st channel length        
#    tlm_set[i].plot_IdVg() #Plotting Id-Vg for all channel lengths in linear scale
#    tlm_set.plot_IdVd(60) #Plotting Id-Vd for all channel lengths @ Vgs=60V
#    tlm_set[i].plot_IdVd(25) #Plotting Id-Vd for all channel lengths @ Vgs=50V
#    tlm_set[i].plot_IdVd(50) #Plotting Id-Vd for all channel lengths @ Vgs=50V
#    print(tlm_set[i].Rc[-1])
#    print(tlm_set[i].del_Rc[-1])
#    print(tlm_set[i].Rc[-1]) #Print minimum contact resistan\ce
#    tlm_set[i].plot_Rsh('n2D') #Plotting Rsh vs Vov    
#    tlm_set[i].plot_Rtot('n2D') #Plotting Rsh vs Vov    
###############################################################################        
    if tlm_set[i].R_squared[-1] > 0.0:#
        
#    if tlm_set[i].R_squared > 0.9:#
#    if tlm_set[i].R_squared>0.856 and tlm_set[i].R_squared<0.857:
#    if file_prefix[i] == '7-4-':
        fig, (ax1, ax2) = plt.subplots(1,2)
#        fig2, ax6 = plt.subplots()
#        tlm_set[i].plot_IdVg(ax1,fig,'log') #Plotting Id-Vg for all channel lengths in log scale
        tlm_set[i].plot_IdVg_Vov(ax1,fig,'log') #Plotting Id-Vg for all channel lengths in log scale
#        tlm_set[i].plot_TLM_fit(ax3,fig,int(tlm_set[i].Vov.shape[0]/2)) #Plotting Rc vs n2D
#        print(tlm_set[i].sbh)
        tlm_set[i].plot_sbh(ax2,fig,0) #Plotting Rsh vs Vov    
#        tlm_set[i].plot_sbh(ax5,fig,1) #Plotting Rsh vs Vov    
#        tlm_set[i].plot_sbh(ax6,fig,2) #Plotting Rsh vs Vov    
#        tlm_set[i].plot_IdVg(ax6,fig2,'log') #Plotting Id-Vg for all channel lengths in log scale
#        tlm_set[i].plot_IdVd(ax2,fig,20) #Plotting Id-Vd for all channel lengths @ Vgs=50V    
#        tlm_set[i].idvd[0].plot_IdVd(ax3,fig) #Plotting Id-Vd for 1st channel length        
#        tlm_set[i].plot_Rc(ax2,fig,'n2D') #Plotting Rc vs n2Ds
#        tlm_set[i].plot_mu_TLM(ax3,fig,'n2D') #Plotting Rc vs n2D            
#        print(tlm_set[i].del_Rc)        
        plt.show()        
#        print(tlm_set[i].R_squared)
#        print(tlm_set[i].Vg_idvg_forward[-2:])
        fig.savefig("summary-"+str(i)+".svg",transparent=True, bbox_inches='tight', pad_inches=0.1)
#        print(str(tlm_set[i].Rc[-1])+" +/- "+str(tlm_set[i].del_Rc[-1]))
    
    if isfolder_organized[i]:
        os.chdir(target_dir)