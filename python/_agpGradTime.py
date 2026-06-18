import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

N = 8
E0 = 1.0
k0 = 1

alphas = np.array([0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4])
alpha = 0.5

fig = plt.figure()

for i in range(len(alphas)):
    alpha = alphas[i]
    FileName = "./../data/alphaInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    df = pd.read_csv(FileName)
    AGPgradNorm = np.array(df["AGPGradNorm"].values.tolist())
    t = np.array(df["Time"].values.tolist())
    plt.plot(t,AGPgradNorm,label=r"$\alpha$ = %.3f"%alpha)

plt.legend()
plt.show()