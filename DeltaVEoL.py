''' Delta V required for end of life manoeuvre '''


import numpy as np
from matplotlib import pyplot as plt


def mass_to_radius(m, rho, t):
    r = np.sqrt(m/9.6/rho/np.pi/t)
    return r


def radius_to_volume(r):
    V = 10/3*np.pi*r**3
    return V


def volume_to_radius(V):
    r = (V*3/10/np.pi)**(1/3)
    return r


def radius_to_mass(r, rho, t):
    m = 9.6*rho*np.pi*r**2*t
    return m

# Constants
first_mass = 100
last_mass = 501
original_dry_masses = np.array([np.arange(first_mass, last_mass), np.arange(first_mass, last_mass), np.arange(first_mass, last_mass)])
dry_masses = np.ones((3, last_mass-first_mass))
for i in range(3):
    for j in range(len(dry_masses[0, :])):
        rango = np.arange(first_mass, last_mass, 1)
        dry_masses[i, j] = rango[j]
initial_orbital_altitude = 1000  # Examples [km]
altitude_for_reentry_dueToDrag = 100  # Example [km]
mass_Earth = 5.972 * (10 ** 24)  # [kg]
Grav_constant = 0.0000000000667  # [m^3/kg/s^2]
g0 = 9.80665  # [m/s^2] gravitational acceleration at Earth's surface
radius_Earth = 6371  # [km]
Prop_density = np.array([1021, 1450, 1036])  # Hydrazine, Hydrogen peroxide, Propyl nitrate
Isp = np.array([283, 117, 210])  # Hydrazine, Hydrogen peroxide, Propyl nitrate
Tank_density = 4430  # Monopropellant tank made of Titanium 6A1-4V
tank_thickness = 1.2*10**(-3)  # [m]
tank_mass0 = 10  # [kg]
r_tank0 = mass_to_radius(tank_mass0, Tank_density, tank_thickness)  # [m]
V_tank0 = radius_to_volume(r_tank0)  # [m^3]

# Iteration parameters
counter = 0
limit = 5
iteration_n = np.arange(1, limit+1, 1)
mass_variation = np.ones(len(iteration_n))

# Delta v computation
tot_M_prop = np.ones((3, len(dry_masses[0, :])))
r = np.array(radius_Earth + initial_orbital_altitude)  # km
a = np.array(0.5 * (radius_Earth + altitude_for_reentry_dueToDrag + r))
orbital_velocity = np.array(np.sqrt((mass_Earth * Grav_constant) / r / 1000))
V_a = np.array(np.sqrt((mass_Earth * Grav_constant) * ((2 / r / 1000) - (1 / a / 1000))))
Delta_V = np.array(orbital_velocity - V_a)

# Iteration

M_prop_previous = np.zeros((3, len(dry_masses[0, :])))

while counter < limit:

    # From delta v to propellant mass using the rocket equation
    M_props = np.ones((3, len(dry_masses[0, :])))
    for j in range(3):
        for idx in range(len(M_props[0, :])):
            M_props[j, idx] = dry_masses[j, idx]*np.exp(Delta_V/Isp[j]/g0) - dry_masses[j, idx]

    # tank calculations
    Delta_volumes = np.ones((3, len(M_props[0, :])))
    for idx in range(len(Prop_density)):
        Delta_volumes[idx, :] = (M_props[idx, :] - M_prop_previous[idx, :])/Prop_density[idx]

    for j in range(3):
        for k in range(len(M_props[0, :])):
            final_volume = Delta_volumes[j, k] + V_tank0
            final_radius = volume_to_radius(final_volume)
            final_mass = radius_to_mass(final_radius, Tank_density, tank_thickness)
            delta_mass = final_mass - tank_mass0  # increase in mass of the tank

            # update dry masses
            dry_masses[j, k] = dry_masses[j, k] + delta_mass
            mass_variation[counter] = dry_masses[j, k]
    counter += 1
    M_prop_previous = M_props

plt.figure()
plt.plot(iteration_n, mass_variation)
plt.xlabel('Iteration')
plt.ylabel('Dry mass [kg]')

plt.figure()
plt.plot(original_dry_masses[0, :], M_props[0, :], original_dry_masses[0, :], dry_masses[0, :])
plt.xlabel('Dry mass')
plt.ylabel('Mass [kg]')
plt.legend(('Propellant mass (Hydrazine)', 'Updated dry mass'))
plt.grid()

plt.figure()
plt.plot(original_dry_masses[1, :], M_props[1, :], original_dry_masses[1, :], dry_masses[1, :])
plt.xlabel('Dry mass')
plt.ylabel('Mass [kg]')
plt.legend(('Propellant mass (Hidrogen Peroxide)', 'Updated dry mass'))
plt.grid()

plt.figure()
plt.plot(original_dry_masses[2, :], M_props[2, :], original_dry_masses[2, :], dry_masses[2, :])
plt.xlabel('Dry mass')
plt.ylabel('Mass [kg]')
plt.legend(('Propellant mass (Propyl Nitrate)', 'Updated dry mass'))
plt.grid()

plt.show()




