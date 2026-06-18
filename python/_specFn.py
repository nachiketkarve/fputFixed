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

#betas = np.array([0.1, 2, 6, 10])
alphas = np.array([0.51])

fig = plt.figure()

for alpha in alphas:
    try:
        FileName = ".//..//data//alphaSpecFn-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
        df = pd.read_csv(FileName)
        f = df['Frequency'].values
        s = df['Spectral Function'].values

        p = np.polyfit(np.log(f[1:int(len(f)/10000)]),np.log(s[1:int(len(f)/10000)]),1)

        plt.plot(f[1:int(len(f))],s[1:int(len(s))])
        plt.plot(f[1:int(len(f)/10000)], np.exp(p[0]*np.log(f[1:int(len(f)/10000)]) + p[1]), 'k',label=r"$S(\omega) \sim \omega^{%.3f}$"%p[0])
        
    except FileNotFoundError:

        continue

plt.xscale('log')
plt.yscale('log')
plt.legend() 

plt.show()