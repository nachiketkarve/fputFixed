# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 20:39:46 2025

@author: karve
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
#sys.path.append("F:\\FPUT\\fput")
#import fputAlpha as fput
import matplotlib.animation as ani
from matplotlib import ticker

N = 8
E0 = 1.0
k0 = 1
alpha = 1

Q1 = 0.2
Q2 = 0.2

FileName = ".//..//alphaInit-N" + str(int(N)) + "-Q1" + "{:.6f}".format(Q1) + "-Q2" + "{:.6f}".format(Q2) + "-A" + "{:.6f}".format(alpha) + ".csv"
df = pd.read_csv(FileName)


Q1s = np.linspace(0,Q1,df.values.shape[0])
Q2s = np.linspace(0,Q2,df.values.shape[1])

Q1Start = 0.04
i1Start = int(Q1Start/(Q1s[1]-Q1s[0]))

Q2Start = 0.0
i2Start = int(Q2Start/(Q2s[1]-Q2s[0]))

plt.figure()
X,Y = np.meshgrid(Q1s[i1Start:len(Q1s)],Q2s[i2Start:len(Q2s)])
cs = plt.pcolormesh(X,Y,np.transpose(np.log10(df.values[i1Start:len(Q1s), i2Start:len(Q2s)])),cmap='plasma')
cbar = plt.colorbar(cs)
plt.xlabel(r"$Q_1$")
plt.ylabel(r"$Q_2$")
cbar.set_label(r"$\sigma^2[\mathcal{A}_\alpha]$")
plt.tight_layout()
#plt.plot(Q1sbreather,Q2sbreather)
plt.savefig("_AGPNorm.pdf",dpi=100)