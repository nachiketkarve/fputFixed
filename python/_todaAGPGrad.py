# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:51:57 2024

@author: karve
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("F:\\FPUT\\python\\toda")
import toda as toda
import matplotlib.animation as ani

N = 4
E0 = 1.0
k0 = 1
alphaMax = 4
alpha = 0.01
deltaAlpha = 0.01
frames = 500

yCircle = np.linspace(-1,1,100)
xCircle = np.sqrt(1-yCircle**2)

A,V,C,FourierComponents = toda.initialize(N)

alphas = np.array([])
gradMag = np.array([])
gradMagX = np.array([])
gradMagP = np.array([])

while alpha <= alphaMax:
    FileName = "./../data/todaBr-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    try:
        df = pd.read_csv(FileName)
        period = df["period"].values.tolist()[0]
        Q = np.array(df["Q"].values.tolist()[0:N+2])
        P = np.zeros(N+2)
        q,p = toda.FT(Q,P,FourierComponents)
        En = toda.TotalEnergy(q,p,alpha)
        AGPgrad = np.array(df["AGPgrad"].values.tolist())/En
        if period < 99.0:
            alphas = np.append(alphas,np.sqrt(En)*alpha)
            gradMag = np.append(gradMag,np.linalg.norm(AGPgrad))
            gradMagX = np.append(gradMagX,np.linalg.norm(AGPgrad[0:N]))
            gradMagP = np.append(gradMagP,np.linalg.norm(AGPgrad[N+1:2*N]))
        
        alpha = alpha + deltaAlpha
    
    except FileNotFoundError:
        alpha = alpha + deltaAlpha


plt.plot(alphas,gradMag)