"""
Created on Mon Apr  6 23:18:02 2020

@author: Aravindh Kumar
"""

"""
Code to analyse cascade/Janis data
"""

# import sys
import os
# Including 3rd-party libraries
# import matplotlib
# matplotlib.rcParams['text.usetex'] = True #For rendering latex in plot elements
import matplotlib.pyplot as plt
import matplotlib as mpl
# from scipy.signal import argrelextrema
from scipy.signal import find_peaks  # For finding local maxima in gm/dgm

# Changing plot properties
plt.style.use(['default', 'seaborn-bright'])
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['xtick.labelsize'] = 9
mpl.rcParams['ytick.labelsize'] = 9
mpl.rcParams['font.size'] = 10
mpl.rcParams['lines.markersize'] = 3
mpl.rcParams['figure.subplot.wspace'] = 0.3  # Horizontal spacing in subplots
mpl.rcParams['figure.subplot.hspace'] = 0.3  # Vertical spacing in subplots
mpl.rcParams['figure.figsize'] = 12, 6  # Figure size
mpl.rcParams['figure.constrained_layout.h_pad'] = 0.04167  # Padding around axes objects
mpl.rcParams['figure.constrained_layout.w_pad'] = 0.04167
mpl.rcParams['legend.fontsize'] = 'small'
# mpl.rcParams['figure.autolayout'] = True

# plt.style.use('default')
import pandas as pd
import xlrd
import numpy as np
import seaborn as sns
# import scipy.linalg
import statsmodels.api as sm
# from statsmodels.api import graphics
from pathlib import Path

a = np.arange(-3,4)
b = np.vstack((-6*a, -4*a, -3*a, -2*a, 2*a, 3*a, 4*a,6*a))
fig, ax = plt.subplots()
ax.plot(a,b.T)
plt.show()