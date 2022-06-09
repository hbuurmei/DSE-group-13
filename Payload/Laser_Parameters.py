import numpy as np

### Constants ###
# Cm0
C_m0 = 420E-6
# Wavelength
Lambda = 532E-9
# Other efficiency losses (atmospheric, etc.)
T_eff = 0.9
# Beam quality factor
M_squared = 2
# Diffraction constant
a = 1.7
# Pulse duration
tau = 1E-10
# Optimum fluence
fluence_opt = 8.5E8 * np.sqrt(tau)
# Coupling coefficient
C_m = C_m0 / (8.5E8 * Lambda)**(1/4)

### Independent variables ###
# Range
z = 250000
# Spot size (conservative, taken as 2 times the 0.0103 arcsec pointing accuracy * range)
d_s = 0.1 + 2 * (z * np.tan((0.0103) * (0.000004848136811095)))
# Velocity difference needed to de-orbit a debris particle
# From 1000 km orbit to 340 km: 174.3957059 m/s
# From 750 km orbit to 340 km: 111.71768422957393 m/s
# From 500 km orbit to 340 km: 44.995583855253244 m/s
delta_V_needed = 174.3957059

### Dependent variables ###
# Effective mirror aperture (diameter)
D_eff = a*M_squared*Lambda*z / (d_s)
# Total mirror aperture (diameter)
D = D_eff / 0.9
# Pulse Energy
E = fluence_opt / (4*(D_eff**2)*T_eff / (np.pi * (M_squared**2)*(a**2)*(Lambda**2)*(z**2)))

### Debris parameters ###
# Time applicable to ablate debris
tt_dv = 50
# Debris area-to-mass ratio
AMR = 0.07946
# Acceleration imparted on the debris
a_deb = delta_V_needed/tt_dv

### More dependent variables ###
# Laser pulse frequency
f = a_deb/(fluence_opt*AMR*C_m)
# Laser output power
output_power = f * E

# Laser single efficiencies
efficiency_conversion = 0.53 * 0.7 * 0.7
# Efficiency of ablation at an incidence angle
maximum_ablation_angle = 20
incidence_angle_efficiency = 1 * np.cos(np.pi/180 * maximum_ablation_angle)**(4/3)

# Input power to the laser
input_power = output_power/(efficiency_conversion*incidence_angle_efficiency)

### Laser budgets ###
m_laser_literature = 2500 #kg [at 125 kW]
m_mirror_literature = 400 #kg [at 3m diameter]

m_laser = input_power/125000 * m_laser_literature
m_mirror = (D/3)**2 * m_mirror_literature
m_total = m_laser + m_mirror

print('Power input: ', input_power, 'W')
print('Pulse energy: ', E, 'J')
print('Pulse frequency: ', f, 'Hz')
print('Pulse duration: ', tau, 's')
print('Laser mass: ', m_laser, 'kg')
print('Mirror mass: ', m_mirror, 'kg')
print('Total mass: ', m_total, 'kg')

