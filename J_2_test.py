import numpy as np
import matplotlib.pyplot as plt

# Constants
R_e = 6378.137e3  # [m]
J_2 = 0.00108263  # [-]
mu = 3.986004418e14  # [m3/s2]


def J_2_RAAN(a, e, i):
    n = np.sqrt(mu / a ** 3)
    RAAN_dot = -1.5*n*R_e**2*J_2*np.cos(np.deg2rad(i))/a**2/(1-e**2)**2
    return RAAN_dot


def J_2_w(a, e, i):
    n = np.sqrt(mu / a ** 3)
    w_dot = 0.75*n*R_e**2*J_2*(4 - 5*(np.sin(np.deg2rad(i)))**2)/a**2/(1-e**2)**2
    return w_dot


e = 0
i = np.arange(0, 180, 1)
a = 1000e3 + R_e

RAAN_dot = J_2_RAAN(a, e, i)*180/np.pi*60*60*24
w_dot = J_2_w(a, e, i)*180/np.pi*60*60*24

plt.plot(i, RAAN_dot)
plt.show()
plt.plot(i, w_dot)
plt.show()
