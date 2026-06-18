# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 14:49:21 2024

@author: karve
"""

#import sys
#sys.path.append('F://FPUT//python//fput')

import numpy as np
import matplotlib.pyplot as plt
#import fputAlpha as br
import pandas as pd
#import toda as br
from scipy.optimize import curve_fit
from matplotlib import ticker

isSingleCol = 1
colFrac = 0.6
htFrac = 0.8

doubleColpt = 246
singleColpt = 510

fig_width_pt = (isSingleCol*singleColpt + (1-isSingleCol)*doubleColpt)
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width = fig_width_pt*inches_per_pt
fig_height = fig_width*golden_mean
fig_size = [fig_width*colFrac,fig_width*colFrac*htFrac]
fig_font = 10
params = {'backend' : 'ps',
          'axes.labelsize': fig_font,
          'font.size': fig_font,
          'legend.fontsize': fig_font,
          'xtick.labelsize': fig_font,
          'ytick.labelsize': fig_font,
          'text.usetex': True,
          'figure.figsize': fig_size,
          'font.family': 'serif',
          'font.serif': 'STIX',
          'mathtext.fontset': 'stix'}

plt.rcParams.update(params)

N = 8
alpha = 1.0

Q1 = 0.2
Q2 = 0.2

FileName = "./../alphaInit-N" + str(int(N)) + "-Q1" + "{:.6f}".format(Q1) + "-Q2" + "{:.6f}".format(Q2) + "-A" + "{:.6f}".format(alpha) + ".csv"
df = pd.read_csv(FileName)

Q1s = np.linspace(0,Q1,df.values.shape[0])
Q2s = np.linspace(0,Q2,df.values.shape[1])

Q1Start = 0.04
i1Start = int(Q1Start/(Q1s[1]-Q1s[0]))

Q2Start = 0.0
i2Start = int(Q2Start/(Q2s[1]-Q2s[0]))

plt.figure()
X,Y = np.meshgrid(Q1s[i1Start:len(Q1s)],Q2s[i2Start:len(Q2s)])
cs = plt.pcolormesh(X,Y,np.transpose(df.values[i1Start:len(Q1s), i2Start:len(Q2s)]),cmap='plasma')
cbar = plt.colorbar(cs)
plt.xlabel(r"$Q_1$")
plt.ylabel(r"$Q_2$")
cbar.set_label(r"$\sigma^2[\mathcal{A}_\alpha]$")
plt.tight_layout()
#plt.plot(Q1sbreather,Q2sbreather)
plt.savefig("_AGPNorm.pdf",dpi=100)



"""
Q1s = np.arange(0.05,0.25,0.01)
Q2s = np.arange(0.01,0.25,0.01)

A,V,C,FourierComponents = br.initialize(N)

norms = np.zeros([len(Q1s),len(Q2s)])

for i in range(len(Q1s)):
    for j in range(len(Q2s)):
        Q1 = Q1s[i]
        Q2 = Q2s[j]

        FileName = "alphaInit-N" + str(int(N)) + "-Q1" + "{:.6f}".format(Q1) + "-Q2" + "{:.6f}".format(Q2) + "-A" + "{:.6f}".format(alpha) + ".csv"
        try:
            df = pd.read_csv(FileName)
            t = df['Time'].values
            dHInt = df['dHInt'].values
            E = np.mean(df['Energy'].values)

            p = np.polyfit(t,dHInt,1)
            slope = p[0]
            #slope = dHInt[len(t)-1]/t[len(t)-1]
            print(slope)

            AGP = slope*t - dHInt
            norms[i,j] = np.std(AGP)**2#/E**3

        except FileNotFoundError:
            continue

Q1sbreather = np.array([])
Q2sbreather = np.array([])

for Q1 in Q1s:
    FileName = "./../alphaBreatherAmplitude/Fixed-N" + str(int(N)) + "-K" + str(int(1)) + "-Q" + "{:.6f}".format(Q1) + "-A" + "{:.6f}".format(alpha) + ".csv"
    try:
        df = pd.read_csv(FileName)
        Q = np.array(df['Q'].values.tolist()[0:N+2])
        print(Q)
        Q1sbreather = np.append(Q1sbreather,Q1)
        Q2sbreather = np.append(Q2sbreather,Q[2])
    except FileNotFoundError:
        continue


print(norms)
print(Q2sbreather)

plt.figure()
X,Y = np.meshgrid(Q1s,Q2s)
cs = plt.pcolormesh(X,Y,np.transpose(norms),cmap='plasma',norm='log')
cbar = plt.colorbar(cs)
plt.xlabel(r"$Q_1$")
plt.ylabel(r"$Q_2$")
cbar.set_label(r"$\sigma^2[\mathcal{A}_\alpha]$")
plt.tight_layout()
#plt.plot(Q1sbreather,Q2sbreather)
plt.savefig("_AGPNorm.pdf",dpi=100)

plt.figure()
plt.plot(Q1sbreather,Q2sbreather)
#plt.yscale("log")
plt.show()
"""

