# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 18:36:04 2025

@author: karve
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
#sys.path.append("F:\\FPUT\\fput")
#import fputAlpha as fput
import matplotlib.animation as ani

N = 8
E0 = 1.0
k0 = 1

betas = np.array([0.1, 2, 6, 10])

fig = plt.figure()

for beta in betas:
    try:
        FileName = ".//..//data//betaCorrInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-B" + "{:.6f}".format(beta) + ".csv"
        df = pd.read_csv(FileName)
        t = df['Time'].values
        corr = df['Correlation'].values

        plt.plot(t[0:int(len(t)/2)],corr[0:int(len(t)/2)])
        
    except FileNotFoundError:
        continue