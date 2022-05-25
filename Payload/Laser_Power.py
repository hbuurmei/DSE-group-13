import numpy as np

# Debris Characteristics

# Debris Mass [kg]
debris_mass = 0.07
# Orbital Height [km]
h = 1000
# Perigee desired [km]
perigee_desired = 200
# Coupling Coefficient [N/W]
C_opt = 30e-6
# mu
mu = 3.986004418e14
# Radius of Earth [km]
R_e = 6371

print('To go from an initial circular orbit of altitude', h, 'km, to an elliptical orbit of perigee', perigee_desired, 'km and apogee', h, 'km:')


# Orbital Velocity @ Original, Circular Orbit
v = np.sqrt(mu/((h+R_e)*1000))

# Semi major axis
a = (perigee_desired + h + 2*R_e)*1000/2

# Orbital velocity at new, elliptical orbit
v_ellipse = np.sqrt(mu*(2/(1000*(R_e + h)) - 1/a))

# Delta v necessary
delta_v = v_ellipse - v
print('The delta v needed is', delta_v, 'm/s')

# Energy necessary
delta_E = 0.5*debris_mass*(v**2 - v_ellipse**2)
print('The kinetic energy delivered to the debris needs to be ', delta_E, 'J')

# Momentum delivered by laser
p = debris_mass * delta_v
delta_E_ablation = - p / C_opt
print('The de-orbit energy, including coupling coefficient, required is', delta_E_ablation, 'J')

# Efficiency
eff = delta_E / delta_E_ablation
print('The coupling efficiency is', eff,)


