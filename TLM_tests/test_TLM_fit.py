# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 02:04:47 2020

@author: aravi
"""

import pandas as pd
import xlrd
import numpy as np
import seaborn as sns
import statsmodels.api as sm
from pathlib import Path
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import find_peaks  # For finding local maxima in gm/dgm
from device_analysis_classes import *

# Changing plot properties
plt.style.use(['default', 'seaborn-bright'])
mpl = plot_properties(mpl)
mpl.rcParams['figure.figsize'] = 4, 3  # Figure size

# Close all plots from previous run(s)
plt.close('all')

annotate_x = 0.05
annotate_y1 = 0.85
annotate_y2 = 0.75
annotate_y22 = 0.725

rtot_label = r'Total Resistance, $R_{tot}$ [k$\Omega.\mu $m]'
lch_label = r'Channel Length,  $L_{ch}$  [$\mu$m]'

def plot_TLM_fit(ax, fig, Rc, del_Rc, R_squared, TLM_fit):
    """Plots TLM fit - regression plot"""            
    plt.gca().set_prop_cycle(None)
    ax.set_xlim(left=0, right=round(L[-1]*1.1/0.2)*0.2)  # Sets lower x-limit to 0
    sns.regplot(x=lch_label, y=rtot_label, data=TLM_fit,\
                    ci=95, truncate=False, ax=ax)         
    # truncate = False plots the regression line to fill x-limits
    start, _ = ax.get_ylim()
    ax.set_ylim(bottom=np.minimum(0, start))  # Sets lower y-limit to 0
    Rc_anno = np.round(Rc, 1)
    delRc_anno = np.round(del_Rc, 1)
    ax.annotate(r'$R_{C}$ = ' + str(Rc_anno) + r' $\pm$ ' + str(delRc_anno) + r' k$\Omega$.$\mu$m',
                xy=(annotate_x, annotate_y1), xycoords='axes fraction')
    # ax.annotate(r'at $n_{2D}$ = ' + str(np.round(n2D[-1] / 1e12, 1)) + r' $\times$ $10^{12} cm^{-2}$',
    #             xy=(annotate_x, annotate_y22), xycoords='axes fraction', fontsize=8.5)        
    ax.annotate(r'R-squared = ' + str(np.round(R_squared, 3)),
                xy=(.6, .1), xycoords='axes fraction', fontsize=8.5)
    # set_ticks(ax,5)
    # set_ticks(ax,5,'x',0.2)


# Setting current directory of file
file_name = 'F:\Google Drive\Python\DeviceAnalysis\TLM_tests\peideye_ITO.csv'
# file_name = 'F:\Google Drive\Python\DeviceAnalysis\TLM_tests\nagoya_1e19.csv'
# file_name = 'F:\Google Drive\Python\DeviceAnalysis\TLM_tests\nagoya_1e20.csv'
##############################################################################

csv_file = open(file_name)
TLM_data = np.asarray(pd.read_csv(file_name)) #Reading TLM data
L = TLM_data[:,0] #Channel lengths (um)
Rtot = TLM_data[:,1]/1e3 #Channel lengths (kohm)
# Rtot = TLM_data[:,2] #Channel lengths (kohm)
# L = np.array([0.27059, 0.369523, 0.45975, 0.560227]) #um
# Rtot = np.array([7.4065, 9.4435, 11.636, 13.8485]) #kohm.um

X = sm.add_constant(L.conj().transpose())  # Adding constant term to the model
model = sm.OLS(Rtot.conj().transpose(), X)

results = model.fit()
Rc = results.params[0] / 2  # Contact resistance (Ohm.um)
del_Rc = results.bse[0] / 2  # Error in Contact resistance (Ohm.um)
Rsh = results.params[1]  # Sheet resistance (Ohm/sq)
del_Rsh = results.bse[1]  # Error in Sheet resistance (Ohm/sq)
R_squared = results.rsquared
d = {rtot_label: Rtot.conj().transpose(),
     lch_label: L.conj().transpose()}  # Creating dictionary of Rtot vs L
TLM_fit = pd.DataFrame(data=d)  # Converting dictionary to Pandas dataframe        

fig, ax = plt.subplots(1, 1)
plot_TLM_fit(ax, fig, Rc, del_Rc, R_squared, TLM_fit)  # Plotting TLM fit for last Vov
plt.show()

fig.savefig("peideye_ITO_analysis.svg", transparent=True, bbox_inches='tight', pad_inches=0.1)
fig.savefig("peideye_ITO_analysis.png", bbox_inches='tight', pad_inches=0.1)