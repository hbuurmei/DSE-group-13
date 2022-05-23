import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


def getPosition(a, e, t, M_0):
    '''
    Find the position of the spacecraft in the Keplerian system.
    @param: a, e, t
    @return: true_anomaly, the true anomaly of the spacecraft at time t
    '''
    n = np.sqrt(mu / a ** 3)  # Mean motion
    M = n*t - M_0
    # Solve the equation numerically
    func = lambda E: E - e * np.sin(E) - M
    init_guess = 0
    E = fsolve(func, init_guess)
    E = E[0]
    # Final equation
    true_anomaly = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))
    return true_anomaly


M_0 = 0
a = 149597870.7e3  # [m]
e = 0.017
mu = 1.3271249e20
t = 0.5*365.25*24*3600
pos_test = getPosition(a, e, t, M_0)

print(np.rad2deg(pos_test))

