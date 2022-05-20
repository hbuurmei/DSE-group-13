import numpy as np
import matplotlib.pyplot as plt

# Define constants
R_e = 6378137  # [m]
I_sp = 300  # [s]
g_0 = 9.80665  # [m/s2]
mu = 3.986004418e14  # [m3/s2]


def slowed_orbit(t_b, h, m_dry):
    # Circular orbit is assumed!
    r = R_e + h  # Radial position
    delV = 2*np.pi*r/t_b  # The delta V w.r.t. original orbit
    V_new = np.sqrt(mu/r) - delV
    t = np.linspace(0, t_b, 100)
    m = m_dry*np.exp(-(mu/r**2-V_new**2/r)/(I_sp*g_0)*(t-t_b))
    m_0 = m_dry*np.exp(t_b*(mu/r**2-V_new**2/r)/(I_sp*g_0))  # Initial mass
    m_p = m_0 - m_dry  # Propellant mass
    return m, t, m_p


t_b = 365.25/2*24*3600
h_range = range(350000, 1000000, 50000)
m_dry_range = range(100, 10000+100, 100)
prop_list = []

for h in h_range:
    for m_dry in m_dry_range:
        _, _, m_p = slowed_orbit(t_b, h, m_dry)
        prop_list.append(m_p)

print('The minimum propellant mass:', '{:.2e}'.format(min(prop_list)), '[kg].')


m_dry = 350  # [kg]
r = R_e + 350000  # [m]
delV = np.linspace(0.1, 100, 1000)  # [m/s]
m_p = m_dry*(np.exp((mu/r**2-(np.sqrt(mu/r)-delV)**2/r)/(I_sp*g_0)*2*np.pi*r/delV)-1)  # [kg]
t_b = 2*np.pi*r/delV
plt.plot(delV, m_p)
plt.show()
plt.plot(t_b, m_p)
plt.show()
