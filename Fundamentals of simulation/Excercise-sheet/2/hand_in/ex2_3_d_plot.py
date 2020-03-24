import numpy as np
import pylab as plt

data_rk2 = np.loadtxt("data_ex3_c.dat", delimiter=", ")
data_rk4 = np.loadtxt("data_ex3_d.dat", delimiter=", ")

time_rk2 = data_rk2[:,0]
energy_rk2 = data_rk2[:,3]
plt.plot(time_rk2, energy_rk2)

time_rk4 = data_rk4[:,0]
energy_rk4 = data_rk4[:,3]
plt.plot(time_rk4, energy_rk4)

plt.xlabel("time")
plt.ylabel("rel. Energy error")
plt.legend(("RK2", "RK4"))
plt.show()