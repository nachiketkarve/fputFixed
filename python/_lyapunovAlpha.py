import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

N = 63
E0 = 6.3
k0 = 1

alphas = np.array([0.25])
alpha = 0.5

fig = plt.figure()

for i in range(len(alphas)):
    alpha = alphas[i]
    FileName = "./../data/alphaLya-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    df = pd.read_csv(FileName)
    t = np.array(df["Time"].values.tolist())
    Lya = np.array(df["MaxLyapunov"].values.tolist())
    t = np.delete(t,0)
    Lya = np.delete(Lya,0)
    plt.plot(t,Lya,label=r"$\alpha$ = %.3f"%alpha)

print(1.0/Lya[-1])

plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.show()