import numpy as np

### Constants ###
C_m0 = 420E-6
Lambda = 532E-9
E_p = 90
T_eff = 0.9
M_squared = 2
a = 1.7

### Desired values ###
z = 250000
d_s = 0.13
tau = 1E-10
delta_V_needed = 214.096274016425

# Dependent variables ###
fluence_opt = 8.5E8 * np.sqrt(tau)
C_m = C_m0 / (8.5E8 * Lambda)**(1/4)
D_eff = a*M_squared*Lambda*z / (d_s)
# fluence_at_debris = 4*E_p*(D_eff**2)*T_eff / (np.pi * (M_squared**2)*(a**2)*(Lambda**2)*(z**2))

E = fluence_opt / (4*(D_eff**2)*T_eff / (np.pi * (M_squared**2)*(a**2)*(Lambda**2)*(z**2)))
print(E)

print('The optimum fluence is', fluence_opt)
#print('The coupling coefficient is', Cm, 'N/W')
print('D_eff is', D_eff, 'm')


### Finding the optimum E and f ###

### Constants ###
power_limit = 20000
efficiency_conversion = 0.37
maximum_ablation_angle = 20
incidence_angle_efficiency = np.cos(np.pi/180 * maximum_ablation_angle)**(4/3)

usable_power = power_limit * efficiency_conversion * incidence_angle_efficiency


f = usable_power / E_p
AMR = 0.076
a_deb = fluence_opt* AMR * C_m * f
tt_dv = delta_V_needed/a_deb

print('The frequency of the laser is', f, 'Hz')
print('The time to impart dV is', tt_dv, 's')
print('The energy in kWh is', power_limit*tt_dv/3600000)

