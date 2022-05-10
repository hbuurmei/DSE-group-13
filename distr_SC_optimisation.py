import numpy as np
import pandas as pd

# Define constants
R_e = 6378137  # [m]
g_0 = 9.80665  # [m/s2]
mu = 3.986004418e14  # [m3/s2]


# Import the reference data
debris_info = pd.read_csv("iridium_cosmos_result.csv")
debris_info = debris_info.loc[debris_info["Name"] == 'Kosmos 2251-Collision-Fragment']
debris_info = debris_info[["Semi-Major-Axis [m]", "Eccentricity", "Argument of periapsis [rad]", "Mean Anomaly [rad]"]]
debris_info = debris_info.to_dict()
debris_info["Removed"] = np.zeros(len(debris_info["Semi-Major-Axis [m]"]))


def getPosition(a, e, t, M_0):
    '''
    Find the position of the spacecraft in the Keplerian system.
    @param: a, e, t
    @return: true_anomaly, the true anomaly of the spacecraft at time t
    '''
    n = np.sqrt(mu / a ** 3)  # Mean motion
    M = n*t - M_0

    # Newton's method to solve the equation
    def eq(E, e, M): return E - e * np.sin(E) - M

    def d_eq(E, e): return 1 - e * np.cos(E)

    err, E_old = 100, M
    while abs(err) > 1e-13:  # Make smaller late if possible
        # print(E_old)
        E_new = E_old - (eq(E_old, e, M) / d_eq(E_old, e))
        err = E_new - E_old
        E_old = E_new
    E = E_new  # Choose final value

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


t = 10*100*60
dt = 1
# Spacecraft variables
a_sc = R_e + (789 + 40)*1000
w_sc = 0
e_sc = 0
M_0_sc = 0
while not np.all(debris_info["Removed"] == 1):
    # Compute spacecraft position
    true_anomaly_sc = getPosition(a_sc, e_sc, t, M_0_sc)
    pos_sc = KeplerToCartesian(a_sc, e_sc, w_sc, true_anomaly_sc)
    # Update space debris position
    for i in range(len(debris_info["Semi-Major-Axis [m]"])):
        if debris_info["Removed"][i] == 0:
            true_anomaly_debris = getPosition(debris_info["Semi-Major-Axis [m]"][i], debris_info["Eccentricity"][i], t, debris_info["Mean Anomaly [rad]"][i])
            pos_debris = KeplerToCartesian(debris_info["Semi-Major-Axis [m]"][i], debris_info["Eccentricity"][i], debris_info["Mean Anomaly [rad]"][i], true_anomaly_debris)
            rel_pos = pos_debris - pos_sc
            abs_distance = np.linalg.norm(rel_pos)
            if abs_distance < 100e3:
                debris_info["Removed"][i] = 1
                print('Time: ', t, 'removed debris fragment')
    t += dt
    print(t/3600)

