import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from matplotlib import pyplot as plt

# Define constants
R_e = 6378.137e3  # [m]
g_0 = 9.80665  # [m/s2]
mu = 3.986004418e14  # [m3/s2]
h_collision = 789e3  # [m]
debris_n = 10
a_collision = R_e + h_collision

# Import the reference data
debris_info = pd.read_csv("iridium_cosmos_result.csv")
debris_info = debris_info.loc[debris_info["Name"] == 'Kosmos 2251-Collision-Fragment']  # Only Kosmos fragments
debris_info = debris_info[["Semi-Major-Axis [m]", "Eccentricity", "Argument of periapsis [rad]", "Mean Anomaly [rad]"]]
debris_info["Removed"] = np.zeros(len(debris_info["Semi-Major-Axis [m]"]))
debris_info = debris_info.loc[debris_info["Semi-Major-Axis [m]"] > a_collision]
index_list = debris_info.index.tolist()
debris_info = debris_info.to_dict()


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
    # Final equation
    true_anomaly = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(E / 2))
    return true_anomaly


def KeplerToCartesian(a, e, w, true_anomaly):
    '''
    Convert the a position in the Keplerian system to a cartesian system
    @param: true_anomaly, the true anomaly of the spacecraft at time t
    '''
    p = a * (1-e**2)
    r = p/(1 + e * np.cos(w))  # radius

    # Compute the Cartesian position vector
    X = r * (np.cos(w + true_anomaly))
    Y = 0
    Z = r * (np.sin(w + true_anomaly))
    return np.array([X, Y, Z]).T


t0 = 20*100*60
t = t0
dt = 10
debris_counter = 0
distance_sc = 110000
# Spacecraft variables
a_sc = R_e + h_collision + distance_sc
w_sc = 0
e_sc = 0
M_0_sc = 0

ts = np.array([])
percentages = np.array([])

#while not np.all(debris_info["Removed"] == 1):
while debris_counter/debris_n < 0.5:
    ts = np.append(ts, t)
    # Compute spacecraft position
    true_anomaly_sc = getPosition(a_sc, e_sc, t, M_0_sc)
    pos_sc = KeplerToCartesian(a_sc, e_sc, w_sc, true_anomaly_sc)
    # Update space debris position
    for i in index_list:
        if debris_info["Removed"][i] == 0:
            true_anomaly_debris = getPosition(debris_info["Semi-Major-Axis [m]"][i], debris_info["Eccentricity"][i], t, debris_info["Mean Anomaly [rad]"][i])
            pos_debris = KeplerToCartesian(debris_info["Semi-Major-Axis [m]"][i], debris_info["Eccentricity"][i], debris_info["Mean Anomaly [rad]"][i], true_anomaly_debris)
            rel_pos = pos_debris - pos_sc
            abs_distance = np.linalg.norm(rel_pos)
            if abs_distance < 100e3:
                debris_info["Removed"][i] = 1
                print('Time (h): ', (t - t0)/3600, 'removed debris fragment')
                debris_counter += 1

    t += dt
    percentages = np.append(percentages, debris_counter/debris_n)

    print("--------------------------------------------")
    print(round((t - t0)/3600, 2))
    if (round(t))%2 == 0:
        print(debris_counter)
        print(str(round(debris_counter/debris_n*100, 2)) + '%')


plt.figure()
plt.plot((ts - t0)/3600, percentages)
plt.xlabel('Time [hr]')
plt.ylabel('Percentage of debris removed [%]')
plt.grid()
plt.show()
