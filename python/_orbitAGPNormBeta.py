# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 14:37:54 2025

@author: karve
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
#sys.path.append("F:\\FPUT\\fput")
#import fputAlpha as fput
import matplotlib as mpl

N = 8
E0 = 1.0
k0 = 1

betas = np.arange(1.8,1.9,0.01)

figName1 = ".//..//plot//_betaInitAGPNorm_N" + str(int(N)) + "_K" + str(int(k0)) + "_E" + "{:.6f}".format(E0) + ".gif"

fig1,ax1 = plt.subplots()
fig2,ax2 = plt.subplots()
fig3,ax3 = plt.subplots()
fig4,ax4 = plt.subplots()

norm = mpl.colors.Normalize(vmin=np.min(E0*betas),vmax=np.max(E0*betas))
cmap = plt.get_cmap('viridis')
s_m = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
s_m.set_array([])

for i in range(len(betas)):
    beta = betas[i]
    FileName = ".//..//data//betaInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-B" + "{:.6f}".format(beta) + ".csv"
    try:
        df = pd.read_csv(FileName)
        t = df["Time"].values
        en = df["Energy"].values
        entropy = df["Entropy"].values
        entropy = entropy * np.log(N)/np.log(N/2.0)
        AGPnorm = df["AGPnorm"].values
        dHInt = df["dHInt"].values
        
        plt.figure(fig1.number)
        plt.plot(t,AGPnorm,color=s_m.to_rgba(beta))
        
        plt.figure(fig2.number)
        plt.plot(t,en,color=s_m.to_rgba(beta))
        
        plt.figure(fig3.number)
        plt.plot(t,entropy,color=s_m.to_rgba(beta))
        
        plt.figure(fig4.number)
        plt.plot(t,dHInt,color=s_m.to_rgba(beta))
    except FileNotFoundError:
        continue


plt.figure(fig1.number)
plt.xlabel(r"$T$")
plt.ylabel(r"$||A_\beta(T)||^2$")
plt.grid()
plt.xlim([np.min(t),np.max(t)])
plt.colorbar(s_m)
plt.tight_layout()

plt.figure(fig2.number)
plt.xlabel(r"$T$")
plt.ylabel(r"$E$")
plt.grid()
plt.xlim([np.min(t),np.max(t)])
plt.colorbar(s_m)
plt.tight_layout()

plt.figure(fig3.number)
plt.xlabel(r"$T$")
plt.ylabel(r"$S$")
plt.grid()
plt.xlim([np.min(t),np.max(t)])
plt.colorbar(s_m)
plt.tight_layout()

plt.figure(fig4.number)
plt.xlabel(r"$T$")
plt.ylabel(r"$\int_0^T d\tau \ \partial_\beta H(\tau)$")
plt.grid()
plt.xlim([np.min(t),np.max(t)])
plt.colorbar(s_m)
plt.tight_layout()