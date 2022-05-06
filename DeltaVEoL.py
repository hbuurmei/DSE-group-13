'''Optimisation of the Delta V required for end of life'''


import numpy as np


##### CONSTANTS #####
initial_orbital_altitude = [350, 400, 450] #Examples [km]
altitude_for_reentry_dueToDrag = 100 #Example [km]
mass_Earth = 5.972 * (10 ** 24) #[kg]
Grav_constant = 0.0000000000667 #[m^3/kg/s^2]
radius_Earth = 6371 #[km]
Prop_density = [1021, 1450, 1036] #Hydrazine, Hydrogen peroxide, Propyl nitrate
Tank_density = 4430 #Monopropellant tank made of Titanium 6A1-4V

### EQUATIONS ###
r = [radius_Earth + initial_orbital_altitude[i] for i in range(len(initial_orbital_altitude))]
a = [0.5 * (radius_Earth + altitude_for_reentry_dueToDrag + r[i]) for i in range(len(r))]
orbital_velocity = [np.sqrt((mass_Earth * Grav_constant) / (initial_orbital_altitude[i] + radius_Earth)) for i in range(len(initial_orbital_altitude))]
V_a = [np.sqrt((mass_Earth * Grav_constant) * ((2 / r[i]) - (1 / a))) for i in range(len(r))]
# Delta_V = [orbital_velocity - V_a [i] for i in range(len(radius_a1))]