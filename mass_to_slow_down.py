import numpy as np
import matplotlib.pyplot as plt

# Define constants
R_e = 6378137  # [m]
I_sp = 300  # [s]
g_0 = 9.80665  # [m/s2]
mu = 3.986004418e14  # [m3/s2]


def slowed_orbit(delV, h, m_dry):
    # Circular orbit is assumed!
    r = R_e + h  # Radial position
    V_new = np.sqrt(mu/r) - delV
    t_b = 2*np.pi*r/delV  # Burn time
    t = np.linspace(0, t_b, 100)
    m = m_dry*np.exp((mu/r**2-V_new**2/2)/(I_sp*g_0)*(t-t_b))
    m_0 = m_dry*np.exp(t_b*(V_new**2/r-mu/r**2)/(I_sp*g_0))  # Initial mass
    m_p = m_0 - m_dry  # Propellant mass
    return m_p, t_b


delV_range = range(10, 1000+10, 10)
h_range = range(350, 1000+50, 50)
m_dry_range = range(100, 10000+100, 100)

m_p_min, t_b_min = slowed_orbit(100, 1000e3, 500)
print(m_p_min, t_b_min)

for delV in delV_range:
    for h in h_range:
        for m_dry in m_dry_range:
            m_p, t_b = slowed_orbit(delV, h, m_dry)
            if m_p > m_p_min:
                m_p_min = m_p
            if t_b < t_b_min:
                t_b_min = t_b

print(m_p_min, t_b_min)
