import numpy as np
import pylab as plt

data = np.loadtxt("data_ex3_c.dat", delimiter=", ")

time = data[:,0]
energy = data[:,3]
plt.plot(time, energy)
plt.show()