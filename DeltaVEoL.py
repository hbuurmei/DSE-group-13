''' Optimisation of the Delta V required for end of life '''


import numpy as np
from matplotlib import pyplot as plt

##### CONSTANTS #####
dry_masses = np.arange(100, 500)
initial_orbital_altitude = np.array([350, 800, 1200]) #Examples [km]
altitude_for_reentry_dueToDrag = 100 #Example [km]
mass_Earth = 5.972 * (10 ** 24) #[kg]
Grav_constant = 0.0000000000667 #[m^3/kg/s^2]
radius_Earth = 6371 #[km]
Prop_density = np.array([1021, 1450, 1036]) #Hydrazine, Hydrogen peroxide, Propyl nitrate
Tank_density = 4430 #Monopropellant tank made of Titanium 6A1-4V
new_dry_mass = 1000
### EQUATIONS ###

# Delta v computation
while new_dry_mass - dry_masses[0] > 1:
    tot_M_prop = np.ones(9, len(dry_masses))
    for altitude_idx in range(len(initial_orbital_altitude)):
        r = np.array(radius_Earth + initial_orbital_altitude[altitude_idx]) #km
        a = np.array(0.5 * (radius_Earth + altitude_for_reentry_dueToDrag + r))
        orbital_velocity = np.array(np.sqrt((mass_Earth * Grav_constant)/r/1000))
        V_a = np.array(np.sqrt((mass_Earth * Grav_constant) * ((2 / r/1000) - (1 / a/1000))))
        Delta_V = np.array(orbital_velocity - V_a)
        print(Delta_V)

        # From delta v to propellant mass
        Isp = np.array([283, 117, 210])
        g0 = 9.80665
        M_props = np.ones((3, len(dry_masses)))
        for idx in range(len(M_props[0, :])):
            M_props[:, idx] = dry_masses[idx]*np.exp(Delta_V/Isp/g0) - dry_masses[idx]
        print(M_props)


        plt.figure(altitude_idx)
        plt.plot(dry_masses, M_props[0, :], dry_masses, M_props[1, :], dry_masses, M_props[2, :])
        plt.legend(('Hydrazine', 'Hidrogen Peroxide', 'Propyl Nitrate'))
        plt.xlabel('Dry mass')
        plt.ylabel('Propellant mass required for reentry')

        Delta_volumes = np.ones((3, len(M_props[0, :])))
        for idx in range(len(Prop_density)):
            Delta_volumes[idx, :] = M_props[idx, :]/Prop_density[idx]

        Delta_dry_mass = Delta_volumes*Tank_density
        dry_masses += Delta_dry_mass
        print(Delta_dry_mass)
plt.show()




