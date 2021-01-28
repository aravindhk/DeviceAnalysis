import os
import matplotlib.pyplot as plt
from scipy.signal import find_peaks  # For finding local maxima in gm/dgm

import pandas as pd
import xlrd
import numpy as np
import seaborn as sns
import statsmodels.api as sm
from pathlib import Path

fig, ax = plt.subplots()
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.plot(0.1,0.5,'o', color='black')
ax.plot(0.1,0.2,'o', color='black')
ax.plot(0.3,0.4,'o', color='black')
ax.plot(0.6,0.4,'o', color='black')

plt.show()