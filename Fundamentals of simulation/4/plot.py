import numpy as np
import pylab as plt

data_comp = np.loadtxt("data_ex4_1.dat", delimiter=", ", skiprows=1)

plt.plot(data_comp[:,0], data_comp[:,1], "x-")
plt.plot(data_comp[:,0], data_comp[:,4], "-")


plt.xlabel("t")
plt.ylabel("x")
plt.title("Comparison of forward Euler with analytical solution")

plt.legend(("Forward Euler", "Analytical"), loc=4)

plt.show()

