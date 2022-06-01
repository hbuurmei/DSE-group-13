import numpy as np

### Constants ###
Cm0 = 420E-6
Lambda = 532E-9
Pulse_Energy = 80
T_eff = 0.9
M_squared = 2
a = 1.7
z = 250000
d_s = 0.10
tau = 8E-9



# Variables
Fluence_opt = 8.5E8 * np.sqrt(tau)
Cm = Cm0 / (8.58E8 * Lambda)**(1/4)
D_eff = a*M_squared*Lambda*z / (d_s)
Fluence = 4*Pulse_Energy*(D_eff**2)*T_eff / (np.pi * (M_squared**2)*(a**2)*(Lambda**2)*(z**2))

print(Fluence_opt)
print(Cm)
# print('The coupling coefficient is', Cm, 'N/W')
print('The aperture diameter at the space craft is', D_eff, 'm')
print('The fluence at the target it', Fluence, 'J/m^2')