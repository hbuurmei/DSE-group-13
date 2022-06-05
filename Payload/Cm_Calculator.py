import numpy as np

### Constants ###
C_m0 = 420E-6
Lambda = 532E-9
E_p = 90
T_eff = 0.9
M_squared = 2
a = 1.7
tau = 1E-10
fluence_opt = 8.5E8 * np.sqrt(tau)
C_m = C_m0 / (8.5E8 * Lambda)**(1/4)

### Independent variables ###
z = 250000
d_s = 0.11
delta_V_needed = 214

# Dependent variables ###
D_eff = a*M_squared*Lambda*z / (d_s)
# fluence_at_debris = 4*E_p*(D_eff**2)*T_eff / (np.pi * (M_squared**2)*(a**2)*(Lambda**2)*(z**2))
print(C_m)

E = fluence_opt / (4*(D_eff**2)*T_eff / (np.pi * (M_squared**2)*(a**2)*(Lambda**2)*(z**2)))
print('Ep =', E)

print('The optimum fluence is', fluence_opt)
#print('The coupling coefficient is', Cm, 'N/W')
print('D_eff is', D_eff, 'm')


### Finding the optimum E and f ###

### Constants ###
tt_dv = 50


AMR = 0.076
a_deb = delta_V_needed/tt_dv
f = a_deb/(fluence_opt*AMR*C_m)
output_power = f * E



efficiency_conversion = 0.53 * 0.7 * 0.7


maximum_ablation_angle = 20
incidence_angle_efficiency = np.cos(np.pi/180 * maximum_ablation_angle)**(4/3)
input_power = output_power/(efficiency_conversion*incidence_angle_efficiency)

# power_limit = 20000
#
# output_power = power_limit * efficiency_conversion * incidence_angle_efficiency
#
# f = output_power / E
# AMR = 0.076
# a_deb = fluence_opt* AMR * C_m * f
# tt_dv = delta_V_needed/a_deb

print('The frequency of the laser is', f, 'Hz')
print('The time to impart dV is', tt_dv, 's')
print('The energy in kWh is', input_power*tt_dv/3600000)

#BUDGETS
m_laser_phipps = 2500 #kg [at 125 kW]
m_mirror_phipps = 400 #kg [at 3m diameter]

m_laser = input_power/125000 * m_laser_phipps
m_mirror = (D_eff/3)**2 * m_mirror_phipps
m_total = m_laser + m_mirror

print('Power input: ', input_power, 'W')
print('Laser mass: ', m_laser, 'kg')
print('Mirror mass: ', m_mirror, 'kg')
print('Total mass: ', m_total, 'kg')

