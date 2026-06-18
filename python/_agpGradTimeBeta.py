import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

N = 8
E0 = 1
k0 = 1

betas = np.array([0.5])

fig = plt.figure()

for i in range(len(betas)):
    beta = betas[i]
    FileName = "./../data/betaInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-B" + "{:.6f}".format(beta) + ".csv"
    df = pd.read_csv(FileName)
    AGPgradNorm = np.array(df["AGPGradNorm"].values.tolist())
    t = np.array(df["Time"].values.tolist())
    plt.plot(t,AGPgradNorm,label=r"$\beta$ = %.3f"%beta)

plt.legend()
plt.show()