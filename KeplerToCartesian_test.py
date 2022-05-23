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


def KeplerToCartesian(a, e, w, true_anomaly, i, RAAN, position):
    '''
    Convert a position in the Keplerian system to a cartesian system
    @param: true_anomaly, the true anomaly of the spacecraft at time t
    '''
    p = a * (1-e**2)
    r = p/(1 + e * np.cos(true_anomaly))  # radius

    # Compute the Cartesian position vector
    X = r*(np.cos(RAAN) * np.cos(w + true_anomaly) - np.sin(RAAN) * np.sin(
        w + true_anomaly) * np.cos(i))
    Y = r*(np.sin(RAAN) * np.cos(w + true_anomaly) + np.cos(RAAN) * np.sin(
        w + true_anomaly) * np.cos(i))
    Z = r*(np.sin(i) * np.sin(w + true_anomaly))
    position[0] = X
    position[1] = Y
    position[2] = Z
    return position


M_0 = 0
a = 149597870.7e3  # [m]
e = 0.017
mu = 1.3271249e20
RAAN = 0
i = 0
w = np.deg2rad(-77.7906)
t = 168*24*3600
true_anomaly = getPosition(a, e, t, M_0)
position = np.zeros([3, 1])
pos_test = KeplerToCartesian(a, e, w, true_anomaly, i, RAAN, position)
print(pos_test)
