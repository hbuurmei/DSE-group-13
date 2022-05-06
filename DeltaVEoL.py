import numpy as np



orbital_altitude = [350000, 400000, 450000] #Examples
altitude_for_reentry_dueToDrag = 100000 #Example
mass_Earth = 5.972 * (10 ** 24)
Grav_constant = 0.0000000000667
radius_Earth = 6371000
radius_p1 = altitude_for_reentry_dueToDrag + radius_Earth
radius_a1 = [radius_Earth + orbital_altitude[i] for i in range(len(orbital_altitude))]

# orbital_velocity = [np.sqrt((mass_Earth * Grav_constant) / (orbital_altitude[i] + radius_Earth)) for i in range(len(orbital_altitude))]
# V_a = [np.sqrt((mass_Earth * Grav_constant) * ((2 / radius_a1[i]) - (1 / a1))) for i in range(len(radius_a1))]
# Delta_V = [orbital_velocity - V_a [i] for i in range(len(radius_a1))]