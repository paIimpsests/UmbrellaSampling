"""
Plotting results obtained from Umbrella Sampling (US) and free energy barrier
calculation scripts
"""
import numpy as np
import matplotlib.pyplot as plt


### Free energy barrier
n = np.empty([])
DGn = np.empty([])
n, DGn = np.loadtxt("DGn_data/DGn_n2000_P17.00.txt", delimiter="\t", unpack=True)

fig, ax =plt.subplots()
ax.scatter(n, DGn, marker="s")
plt.xlabel(r"$\textrm{Cluster size}~n$")
plt.ylabel(r"$\beta \Delta G(n)$")
plt.show()



### Max cluster size evolution
i = np.empty([])
n_i = np.empty([])
i, n_i = np.loadtxt("n_data/n_n2000_P17.00.txt", delimiter="\t", unpack=True)

fig, ax =plt.subplots()
#ax.plot(i, n_i, marker="o", linestyle="")
ax.bar(i, n_i)
plt.xlabel(r"$\textrm{trajectory iteration}~i$")
plt.ylabel(r"$\textrm{Max cluster size}~n$")
plt.show()
