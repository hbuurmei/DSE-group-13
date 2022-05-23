import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from matplotlib import pyplot as plt

# Define constants
R_e = 6378.137e3  # [m]
g_0 = 9.80665  # [m/s2]
J_2 = 0.00108263  # [-]
mu = 3.986004418e14  # [m3/s2]
h_collision = 789e3  # [m]
debris_n = 1000  # total n is 21627 fragments, change this number for simulation speed
a_collision = R_e + h_collision

# Import the reference data
debris_info = pd.read_csv("iridium_cosmos_result.csv")
debris_info = debris_info.loc[debris_info["Name"] == 'Kosmos 2251-Collision-Fragment']  # Only Kosmos fragments
debris_info = debris_info[["Semi-Major-Axis [m]", "Eccentricity", "Inclination [rad]",
                           "Longitude of the ascending node [rad]", "Argument of periapsis [rad]", "Mean Anomaly [rad]"]]
debris_info["Removed"] = np.zeros(len(debris_info["Semi-Major-Axis [m]"]))
debris_info = debris_info.loc[debris_info["Semi-Major-Axis [m]"] > a_collision - 60e3]
debris_info = debris_info.head(debris_n)
index_list = debris_info.index.tolist()
debris_info = debris_info.to_numpy()


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
    X = r * (np.cos(RAAN) * np.cos(w + true_anomaly) - np.sin(RAAN) * np.sin(
        w + true_anomaly) * np.cos(i))
    Y = r * (np.sin(RAAN) * np.cos(w + true_anomaly) + np.cos(RAAN) * np.sin(
        w + true_anomaly) * np.cos(i))
    Z = r * (np.sin(i) * np.sin(w + true_anomaly))
    position[0] = X
    position[1] = Y
    position[2] = Z
    return position


def getVelocity(a, e, w, true_anomaly, i, RAAN, position):
    p = a * (1-e**2)
    r = p/(1 + e * np.cos(true_anomaly))  # radius
    h = np.sqrt(mu*p)
    V_X = (position[0] * h * e / (r * p)) * np.sin(true_anomaly) - (h / r) * (np.cos(RAAN) * np.sin(w + true_anomaly) + np.sin(RAAN) * np.cos(w + true_anomaly) * np.cos(i))
    V_Y = (position[1] * h * e / (r * p)) * np.sin(true_anomaly) - (h / r) * (np.sin(RAAN) * np.sin(w + true_anomaly) - np.cos(RAAN) * np.cos(w + true_anomaly) * np.cos(i))
    V_Z = (position[2] * h * e / (r * p)) * np.sin(true_anomaly) + (h / r) * (np.cos(w + true_anomaly) * np.sin(i))
    V = np.array([V_X, V_Y, V_Z])
    return V


def J_2_RAAN(a, e, i):
    n = np.sqrt(mu / a ** 3)
    RAAN_dot = -1.5*n*R_e**2*J_2*np.cos(i)/a**2/(1-e**2)**2
    return RAAN_dot


def J_2_w(a, e, i):
    n = np.sqrt(mu / a ** 3)
    w_dot = 0.75*n*R_e**2*J_2*(4 - 5*(np.sin(i))**2)/a**2/(1-e**2)**2
    return w_dot


t0 = 72*100*60
t = t0
dt = 50
debris_counter = 0
distance_sc = 40e3
# Spacecraft variables
a_sc = R_e + h_collision + distance_sc
w_sc = 0
e_sc = 0
i_sc = np.average(debris_info[:, 2])
RAAN_sc = np.average(debris_info[:, 3])
M_0_sc = 0

ts = np.array([])
percentages = np.array([])
position_sc = np.zeros([3, 1])
position_debris = np.zeros([3, 1])

# J_2 effect sc
RAAN_dot_sc = J_2_RAAN(a_sc, e_sc, i_sc)
RAAN_drift_sc = RAAN_dot_sc*dt

w_dot_sc = J_2_w(a_sc, e_sc, i_sc)
w_drift_sc = w_dot_sc*dt

# J_2 effect debris
RAAN_dot = J_2_RAAN(debris_info[:, 0], debris_info[:, 1], debris_info[:, 2])
RAAN_drift = RAAN_dot*dt

w_dot = J_2_w(debris_info[:, 0], debris_info[:, 1], debris_info[:, 2])
w_drift = w_dot*dt

while debris_counter/debris_n < 0.822231:
    ts = np.append(ts, t)
    # Update RAAN and w due to J_2 (sc)
    RAAN_sc += RAAN_drift_sc
    w_sc += w_drift_sc
    # Update RAAN and w due to J_2 (debris)
    debris_info[:, 3] += RAAN_drift
    debris_info[:, 4] += w_drift
    # Compute spacecraft position
    true_anomaly_sc = getPosition(a_sc, e_sc, t, M_0_sc)
    pos_sc = KeplerToCartesian(a_sc, e_sc, w_sc, true_anomaly_sc, i_sc, RAAN_sc, position_sc)
    # Update space debris position
    for i in range(len(debris_info[:, 0])):
        if debris_info[i, 6] == 0:
            true_anomaly_debris = getPosition(debris_info[i, 0], debris_info[i, 1], t, debris_info[i, 5])
            pos_debris = KeplerToCartesian(debris_info[i, 0], debris_info[i, 1], debris_info[i, 4], true_anomaly_debris,
                                           debris_info[i, 2], debris_info[i, 3], position_debris)
            rel_pos = pos_debris - pos_sc
            abs_distance = np.linalg.norm(rel_pos)
            if abs_distance < 100e3:
                vel_sc = getVelocity(a_sc, e_sc, w_sc, true_anomaly_sc, i_sc, RAAN_sc, position_sc)
                pointing_laser = - vel_sc
                if sum(pointing_laser*rel_pos) / (np.linalg.norm(pointing_laser) * np.linalg.norm(rel_pos)) > 0.5:
                    debris_info[i, 6] = 1
                    debris_counter = len(debris_info[debris_info[:, 6] == 1])
                else:
                    print("Ineffective geometry")

    t += dt
    percentages = np.append(percentages, debris_counter/debris_n)

    print("--------------------------------------------")
    print(round((t - t0)/3600, 2))
    if (round(t))%2 == 0:
        print(debris_counter)
        print(str(round(debris_counter/debris_n*100, 2)) + '%')


plt.figure()
plt.plot((ts - t0)/3600/24, percentages*100*0.6)
plt.xlabel('Time since start of mission [days]')
plt.ylabel('Percentage of debris removed [%]')
plt.grid()
plt.savefig('Debris_removed.png')
plt.show()

