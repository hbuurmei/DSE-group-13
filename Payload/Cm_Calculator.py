import numpy as np

### Constants ###
C_m0 = 420E-6
Lambda = 532E-9
E_p = 80
T_eff = 0.9
M_squared = 2
a = 1.7

### Desired values ###
z = 250000
d_s = 0.11
tau = 8E-9
delta_V_needed = 243

# Dependent variables ###
fluence_opt = 8.5E8 * np.sqrt(tau)
C_m = C_m0 / (8.58E8 * Lambda)**(1/4)
D_eff = a*M_squared*Lambda*z / (d_s)
fluence_at_debris = 4*E_p*(D_eff**2)*T_eff / (np.pi * (M_squared**2)*(a**2)*(Lambda**2)*(z**2))

print('The optimum fluence is', fluence_opt)
#print('The coupling coefficient is', Cm, 'N/W')
print('D_eff is', D_eff, 'm')
print('The fluence at the target it', fluence_at_debris, 'J/m^2')

### Finding the optimum E and f ###

### Constants ###
power_limit = 29364.13333
efficiency_conversion = 0.37
maximum_ablation_angle = 20
incidence_angle_efficiency = np.cos(np.pi/180 * maximum_ablation_angle)**(4/3)


usable_power = power_limit * efficiency_conversion * incidence_angle_efficiency


f = usable_power / E_p
AMR = 0.08
a_deb = fluence_at_debris * AMR * C_m * f
tt_dv = delta_V_needed/a_deb
print('The time to impart dV is', tt_dv, 's')

