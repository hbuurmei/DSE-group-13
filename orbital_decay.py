import numpy as np
import USSA76
import matplotlib.pyplot as plt
import scipy.integrate as spi
import scipy.optimize as opt

GM = 3.986004418 * 10 ** 5  # km^3/s^2
R_earth = 6371  # km
Cd_sphere = 2.667
seconds_in_year = 365*24*3600


def EoM(t, y, Cd, area_to_mass):
    alt = np.linalg.norm(y[0:2])
    v = np.linalg.norm(y[2:])
    rho = USSA76.calculate_density(alt - R_earth)
    dr_dt = y[2:]
    dv_dt = -GM * y[0:2]/alt**3 - 0.5 * rho * Cd * area_to_mass * y[2:] * v * 1000
    return np.concatenate((dr_dt, dv_dt))


def calculate_final_orbit_accurate(area_to_mass, t_final, alt_0, Cd):
    return spi.solve_ivp(fun=EoM, args=(Cd, area_to_mass), y0=np.array([alt_0, 0, 0, np.sqrt(GM/alt_0)]),
                         t_span=(0, t_final), method="BDF", atol=10e-8, rtol=10e-8)


def calculate_deorbiting_time_approximation(area_to_mass, alt_0, Cd):
    t = 0
    index = np.argmin(USSA76.h <= alt_0) - 1

    rho = USSA76.calculate_density(alt_0 - R_earth)
    delta_a = - 2 * np.pi * Cd * area_to_mass * rho * alt_0 ** 2 * 1000
    delta_L = - USSA76.H[index] / delta_a
    T_orbit = 2 * np.pi * alt_0 ** 1.5 / np.sqrt(GM)
    delta_T = delta_L * T_orbit
    t += delta_T

    for i in range(0, index):
        rho = USSA76.calculate_density(USSA76.h[i + 1] - R_earth)
        delta_a = - 2 * np.pi * Cd * area_to_mass * rho * USSA76.h[i + 1]**2 * 1000
        delta_L = - USSA76.H[i] / delta_a
        T_orbit = 2 * np.pi * USSA76.h[i + 1]**1.5 / np.sqrt(GM)
        delta_T = delta_L * T_orbit
        t += delta_T

    return t/seconds_in_year


def equation_am(area_to_mass, t_final, alt_0, Cd):
    t = calculate_deorbiting_time_approximation(area_to_mass=area_to_mass, alt_0=alt_0, Cd=Cd)
    return t - t_final


def find_area_to_mass(alt_0, Cd, deorbiting_time, x0):
    return opt.fsolve(equation_am, x0=np.array([x0]), args=(deorbiting_time, alt_0, Cd), xtol=1e-10,
                      full_output=True)


# deorbiting_time = 5
# area_to_mass = []
# alt_range = np.linspace(300, 1000, 1000)
# for alt in alt_range:
#     area_to_mass.append(find_area_to_mass(alt_0=alt + R_earth, Cd=Cd_sphere, deorbiting_time=deorbiting_time,
#                                           x0=0.0001))
#
# plt.plot(alt_range, area_to_mass)
# plt.xlabel(f"Initial altitude (km)")
# plt.ylabel(f"Area to mass ratio (m^2/kg) required for deorbiting time of {deorbiting_time} years")
# plt.show()


expanded_rho = 1
deorbiting_time = 1
foam_mass = []
alt_range = np.linspace(300, 1000, 50)
particles_masses = np.genfromtxt('iridium_cosmos_result.csv', delimiter=',', skip_header=1, usecols=6)
sattelite = np.genfromtxt('iridium_cosmos_result.csv', delimiter=',', skip_header=1, usecols=1, dtype='str')
inclination = np.genfromtxt('iridium_cosmos_result.csv', delimiter=',', skip_header=1, usecols=-4)
# a = np.genfromtxt('iridium_cosmos_result.csv', delimiter=',', skip_header=1, usecols=-6)


particles_masses = particles_masses[sattelite == "Kosmos 2251-Collision-Fragment"]
inclination = inclination[sattelite=="Kosmos 2251-Collision-Fragment"]
# a = a[sattelite=="Kosmos 2251-Collision-Fragment"]

particles_masses = np.sort(particles_masses, kind='mergesort')
# inclination = np.sort(inclination, kind='mergesort')
# a = np.sort(a, kind='mergesort')


particles_masses = particles_masses[int(particles_masses.size/4):int(particles_masses.size * 0.75)]
inclination = inclination[int(inclination.size/2):]
print(np.amax(inclination)*180/np.pi, np.min(inclination)*180/np.pi)
# a = a[int(a.size/2):]


def equation_r(r, a_to_m, rho, m_debris):
    return (np.pi * r ** 2) / (4 / 3 * np.pi * r ** 3 * rho + m_debris) - a_to_m


def find_radius(a_to_m, rho, m_debris, x0):
    return opt.fsolve(equation_r, x0=np.array([x0]), args=(a_to_m, rho, m_debris), xtol=1e-10,
                      full_output=True)


# for alt in alt_range:
#     sol = find_area_to_mass(alt_0=alt + R_earth, Cd=Cd_sphere, deorbiting_time=deorbiting_time, x0=0.01)
#     # print(f"Solving for A/m for alt_0 = {alt}: {sol}")
#     area_to_mass = sol[0]
#     m_total = 0
#     for mass in particles_masses:
#         sol = find_radius(a_to_m=area_to_mass, rho=expanded_rho, m_debris=mass, x0=0.0001)
#         # print(f"Solving for r_foam for m  = {mass}, alt_0 = {alt}: {sol}")
#         r = sol[0]
#         m_foam = 4/3 * np.pi * r**3 * expanded_rho
#         m_total += m_foam
#     foam_mass.append(m_total)


# plt.plot(alt_range, foam_mass)
# plt.xlabel(f"Initial altitude (km)")
# plt.ylabel(f"Foam mass required for deorbiting time of {deorbiting_time} year(s)")
# plt.show()

a_to_m = []
for alt in alt_range:
    sol = find_area_to_mass(alt_0=alt + R_earth, Cd=Cd_sphere, deorbiting_time=deorbiting_time, x0=0.01)
    # print(f"Solving for A/m for alt_0 = {alt}: {sol}")
    a_to_m.append(sol[0])
plt.plot(alt_range, a_to_m)
plt.xlabel("Initial altitude [km]")
plt.ylabel(f"Area-to-mass [m^2/kg] ratio required for deorbiting time of {deorbiting_time} year(s)")
plt.show()

