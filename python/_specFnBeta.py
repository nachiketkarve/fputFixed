import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
#sys.path.append("F:\\FPUT\\fput")
#import fputAlpha as fput
import matplotlib.animation as ani

N = 64
E0 = 1.0
k0 = 1

beta = 5.0

fig = plt.figure()


FileName = ".//..//data//betaSpecFn-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-B" + "{:.6f}".format(beta) + ".csv"
df = pd.read_csv(FileName)
f = df['Frequency'].values
s = df['Spectral Function'].values

p = np.polyfit(np.log(f[1:int(len(f)/100)]),np.log(s[1:int(len(f)/100)]),1)

plt.plot(f[1:int(len(f))],s[1:int(len(s))])
plt.plot(f[1:int(len(f)/100)], np.exp(p[0]*np.log(f[1:int(len(f)/100)]) + p[1]), 'k',label=r"$\omega^{%.3f}$"%p[0])
        

plt.xscale('log')
plt.yscale('log')
plt.legend() 

plt.show()