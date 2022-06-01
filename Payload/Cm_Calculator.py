import numpy as np

### Constants ###
Cm0 = 420E-6
Lambda = 532E-9
E_p = 80
T_eff = 0.9
M_squared = 2
a = 1.7
z = 250000
d_s = 0.11
tau = 8E-9
delta_V_needed = 243



# Variables
Fluence_opt = 8.5E8 * np.sqrt(tau)
Cm = Cm0 / (8.58E8 * Lambda)**(1/4)
D_eff = a*M_squared*Lambda*z / (d_s)
Fluence_at_debris = 4*E_p*(D_eff**2)*T_eff / (np.pi * (M_squared**2)*(a**2)*(Lambda**2)*(z**2))

# print('The optimum fluence you would like is', Fluence_opt)
# print('The coupling coefficient is', Cm, 'N/W')
# print('The aperture diameter at the space craft is', D_eff, 'm')
# print('The fluence at the target it', Fluence_at_debris, 'J/m^2')

# Finding the optimum E and f
Power_limit = 29364.13333
Usable_Power = Power_limit*np.cos(np.pi/180 * 20)**(4/3)*0.37
print(Usable_Power)
f = Usable_Power / E_p
AMR = 0.08
a_deb = Fluence_at_debris * AMR * Cm * f
time_to_ablate = delta_V_needed/a_deb
print(time_to_ablate)

