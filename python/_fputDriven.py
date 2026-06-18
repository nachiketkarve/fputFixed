import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

N = 64
E0 = 1.0
beta = 5.0
lmd = 0.1
k0 = 1

Ts = np.arange(1000, 10000, 1000)

TsPlot = np.array([])
deltaEs = np.array([])
vs = np.array([])

for T in Ts:
    try:
        
        FileName = "./../data/betaDriven-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-B" + "{:.6f}".format(beta) + "-T" + "{:.6f}".format(T) + "-L" + "{:.6f}".format(lmd) + ".csv"
        df = pd.read_csv(FileName)
        Ei = np.array(df["Ei"].values.tolist())
        Ef = np.array(df["Ef"].values.tolist())
        dE = np.array(df["deltaE"].values.tolist())
        TsPlot = np.append(TsPlot, T)
        deltaEs = np.append(deltaEs, np.std(dE)**2)
        vs = np.append(vs, lmd / T)
    except:
        continue

print(deltaEs)

mus = 1/TsPlot
p = np.polyfit(np.log(mus), np.log(deltaEs/vs**2), 1)


plt.figure()
plt.plot(1/TsPlot, deltaEs/vs**2)
plt.plot(1/TsPlot, np.exp(p[1]) * mus**p[0], label=r"$\sim \mu^{%.3f}$"%p[0])
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.show()