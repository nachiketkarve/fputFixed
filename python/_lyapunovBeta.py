import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

N = 32
E0 = 1.0
k0 = 1

betas = np.array([0.6,0.9])
beta = 0.5

fig = plt.figure()

for i in range(len(betas)):
    beta = betas[i]
    FileName = "./../data/betaLya-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-B" + "{:.6f}".format(beta) + ".csv"
    df = pd.read_csv(FileName)
    t = np.array(df["Time"].values.tolist())
    Lya = np.array(df["MaxLyapunov"].values.tolist())
    t = np.delete(t,0)
    Lya = np.delete(Lya,0)
    plt.plot(t,Lya,label=r"$\beta$ = %.3f"%beta)
    print(1.0/Lya[-1])

plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()